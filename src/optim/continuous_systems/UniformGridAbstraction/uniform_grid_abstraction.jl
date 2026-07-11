export UniformGridAbstraction

module UniformGridAbstraction

import Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const SY = Dionysos.Symbolic
const OP = Dionysos.Optim
const OPDS = OP.DiscreteSystems

import StaticArrays: SVector, SMatrix
import MathematicalSystems as MS
import HybridSystems
import LazySets
using JuMP
import LinearAlgebra

import Distributed

export Optimizer

include("alternating_simulation_problem.jl")
include("lifted_control_optimizer.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")
include("reach_and_stay_problem.jl")
include("cosafe_ltl_problem.jl")

"""
    Optimizer{T} <: Dionysos.Optim.CompositeDionysosOptimizer

A high-level abstraction-based solver that automatically orchestrates system abstraction and control synthesis.
This wrapper follows the **classical abstraction pipeline** (e.g., as in SCOTS), where the state and input spaces are discretized into hyper-rectangular cells, independent of the specific control task.

It delegates responsibility to modular sub-solvers: one for abstraction and one for control, depending on the type of problem to be solved.

---

### Structure and Sub-solvers

The optimizer internally manages two sub-solvers:

- `abstraction_solver`: [`OptimizerAlternatingSimulationProblem`](@ref):  
  Used to compute the symbolic abstraction of the system from its dynamics and domain.

- `control_solver`: One of the following control-specific optimizers, depending on the problem type:
    - [`OptimizerOptimalControlProblem`](@ref): for [`reachability/optimal control problems`](@ref Dionysos.Problem.OptimalControlProblem).
    - [`OptimizerSafetyProblem`](@ref): for [`safety problems`](@ref Dionysos.Problem.SafetyProblem).

---

### Behavior

- The user sets the control task via:  
  `MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), my_problem)`  
  where `my_problem` is a subtype of [`ProblemType`](@ref Dionysos.Problem.ProblemType).

- The optimizer automatically dispatches to the appropriate control solver based on the problem type.

- If the abstraction has not yet been computed, it is automatically built **before** solving the control problem.

- Once computed, the abstraction is cached — switching the control problem (e.g., from safety to reachability) does not recompute it.

- The field `solve_time_sec` tracks the runtime of the last call to `MOI.optimize!`.

- The resulting controller and value function are stored and can be queried from the wrapper.

---

### User-settable and access to subsolver fields

Via `MOI.set(...)`, the user may configure `abstraction_solver` and `control_solver` parameters.
Any field accessible in the sub-solvers (abstraction or control) can be transparently accessed via:

```julia
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), grid)
MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"))
```

### Example

```julia
using Dionysos, JuMP
optimizer = MOI.instantiate(Dionysos.Optim.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), my_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer)

time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
value_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"))
controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
```
"""
mutable struct Optimizer{T} <: OP.AbstractionControlOptimizer
    abstraction_solver::Union{Nothing, OptimizerAlternatingSimulationProblem{T}}
    control_solver::Union{Nothing, MOI.AbstractOptimizer}
    concrete_controller::Union{Nothing, ST.AbstractContinuousController}
    out_of_domain_handler::SY.AbstractOutOfDomainHandler
    solve_time_sec::T
    print_level::Int

    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, SY.NoOutOfDomainHandler(), 0.0, 1)
    end
end
Optimizer() = Optimizer{Float64}()

OP.default_abstraction_solver(::Optimizer) = OptimizerAlternatingSimulationProblem()

"""
    control_solver_for(problem)

Return a fresh control sub-solver matching the type of the concrete `problem`. Supporting a new
`ProblemType` in this solver family means adding one method here.
"""
control_solver_for(::PR.OptimalControlProblem) = OptimizerOptimalControlProblem()
control_solver_for(::PR.SafetyProblem) = OptimizerSafetyProblem()
control_solver_for(::PR.ReachAndStayProblem) = OptimizerReachAndStayProblem()
control_solver_for(::PR.CoSafeLTLProblem) = OptimizerCoSafeLTLProblem()
control_solver_for(problem) = error("Unsupported problem type: $(typeof(problem))")

function OP.set_concrete_problem!(
    model::Optimizer,
    problem::PR.AlternatingSimulationProblem,
)
    old_abstraction_solver = model.abstraction_solver
    model.abstraction_solver = OptimizerAlternatingSimulationProblem()
    if old_abstraction_solver !== nothing
        model.abstraction_solver.execution_backend =
            old_abstraction_solver.execution_backend
        model.abstraction_solver.print_level = old_abstraction_solver.print_level
        model.abstraction_solver.progress_update_interval =
            old_abstraction_solver.progress_update_interval
        model.abstraction_solver.progress_dt = old_abstraction_solver.progress_dt
    end
    model.abstraction_solver.alternating_simulation_problem = problem
    model.control_solver = nothing
    return
end

function OP.set_concrete_problem!(model::Optimizer, problem)
    model.control_solver = control_solver_for(problem)
    model.control_solver.concrete_problem = problem
    if model.abstraction_solver.alternating_simulation_problem === nothing
        model.abstraction_solver.alternating_simulation_problem =
            PR.AlternatingSimulationProblem(problem.system, nothing)
    end
    return
end

function reset!(optimizer::Optimizer)
    optimizer.concrete_controller = nothing
    optimizer.solve_time_sec = 0.0
    if optimizer.control_solver !== nothing
        reset!(optimizer.control_solver)
    end
    if optimizer.abstraction_solver !== nothing
        reset!(optimizer.abstraction_solver)
    end
    return optimizer
end

function OP.build_concrete_controller(
    model::Optimizer,
    abstract_system,
    abstract_controller,
)
    return SY.quantize_controller(
        SY.without_metadata(abstract_system), # light version of the abstract system, without metadata, for faster controller saving/loading
        abstract_controller;
        out_of_domain_handler = model.out_of_domain_handler,
    )
end

end # module
