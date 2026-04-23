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
using JuMP
import LinearAlgebra

import Distributed

export Optimizer

include("alternating_simulation_problem.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")
include("cosafe_ltl_problem.jl")

"""
    Optimizer{T} <: MOI.AbstractOptimizer

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
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    abstraction_solver::Union{Nothing, OptimizerAlternatingSimulationProblem{T}}
    control_solver::Union{Nothing, MOI.AbstractOptimizer}
    concrete_controller::Union{Nothing, ST.AbstractContinuousController}
    solve_time_sec::T
    print_level::Int

    function Optimizer{T}() where {T}
        return new{T}(nothing, nothing, nothing, 0.0, 1)
    end
end
Optimizer() = Optimizer{Float64}()

MOI.is_empty(optimizer::Optimizer) = optimizer.abstraction_solver === nothing

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    return model.print_level = value ? 0 : 1
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    param_symbol = Symbol(param.name)

    if model.abstraction_solver === nothing
        model.abstraction_solver = OptimizerAlternatingSimulationProblem()
    end

    if param_symbol == :abstract_system
        MOI.set(
            model.abstraction_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            value,
        )
        return
    end

    if param_symbol == :concrete_problem
        # Assign appropriate control solver
        if isa(value, Dionysos.Problem.AlternatingSimulationProblem)
            model.abstraction_solver = OptimizerAlternatingSimulationProblem()
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("alternating_simulation_problem"),
                value,
            )
            model.control_solver = nothing  # No control solver
        elseif isa(value, Dionysos.Problem.OptimalControlProblem)
            model.control_solver = OptimizerOptimalControlProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        elseif isa(value, Dionysos.Problem.SafetyProblem)
            model.control_solver = OptimizerSafetyProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        elseif isa(value, Dionysos.Problem.CoSafeLTLProblem)
            model.control_solver = OptimizerCoSafeLTLProblem()
            MOI.set(
                model.control_solver,
                MOI.RawOptimizerAttribute("concrete_problem"),
                value,
            )
        else
            error("Unsupported problem type: $(typeof(value))")
        end

        # Instantiate an abstraction_solver if it has not already been created
        if model.abstraction_solver.alternating_simulation_problem === nothing
            alternating_simulation_problem =
                Dionysos.Problem.AlternatingSimulationProblem(value.system, nothing)
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("alternating_simulation_problem"),
                alternating_simulation_problem,
            )
        end
        return
    end

    # If the attribute exists in the main optimizer, set it there
    if hasfield(typeof(model), param_symbol)
        return setproperty!(model, param_symbol, value)
    end

    # Try setting it in the sub-solvers
    for solver in (model.abstraction_solver, model.control_solver)
        if solver !== nothing && hasfield(typeof(solver), param_symbol)
            return setproperty!(solver, param_symbol, value)
        end
    end

    return error(
        "Attribute $(param.name) is not recognized by the solver for the considered control problem",
    )
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    param_symbol = Symbol(param.name)

    # First, check if the attribute exists directly in the main optimizer
    if hasfield(typeof(model), param_symbol)
        return getproperty(model, param_symbol)
    end

    # If not found, try to get it from the abstraction solver
    if model.abstraction_solver !== nothing &&
       hasfield(typeof(model.abstraction_solver), param_symbol)
        return getproperty(model.abstraction_solver, param_symbol)
    end

    # If not found, try to get it from the control solver
    if model.control_solver !== nothing &&
       hasfield(typeof(model.control_solver), param_symbol)
        return getproperty(model.control_solver, param_symbol)
    end

    # If the attribute is not recognized, raise an error
    return error(
        "Attribute $(param.name) is not recognized by the solver for the considered control problem",
    )
end

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time_sec
end

function is_abstraction_computed(optimizer::Optimizer)
    return optimizer.abstraction_solver !== nothing &&
           optimizer.abstraction_solver.abstract_system !== nothing
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

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    optimizer.abstraction_solver === nothing &&
        error("The concrete problem is not defined.")

    if !is_abstraction_computed(optimizer)
        MOI.set(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("print_level"),
            optimizer.print_level,
        )
        MOI.optimize!(optimizer.abstraction_solver)
    end

    if optimizer.control_solver !== nothing
        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("print_level"),
            optimizer.print_level,
        )

        abstract_system = MOI.get(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
        )

        MOI.set(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_system"),
            abstract_system,
        )

        MOI.optimize!(optimizer.control_solver)

        abstract_controller = MOI.get(
            optimizer.control_solver,
            MOI.RawOptimizerAttribute("abstract_controller"),
        )

        optimizer.concrete_controller =
            SY.quantize_controller(abstract_system, abstract_controller)
    end

    optimizer.solve_time_sec = time() - t_ref
    return
end

"""
make_out_of_domain_handler(; mode=0, warn=true, dims=nothing)

mode = 0: return nothing when x is not allowed (or outside mapping)
mode = 1: project to nearest allowed abstract state (using mapping coords)
"""
function make_out_of_domain_handler(; mode::Int = 0, warn::Bool = true, dims = nothing)
    if mode == 0
        return (x, abs_sys) -> begin
            Xmap = SY.get_state_mapping(abs_sys)
            Xset = SY.get_state_domain(abs_sys)
            q = MP.get_state_by_coord(Xmap, x)
            if q === nothing || !MP.contains_state(Xset, Xmap, q)
                warn && @warn("State out of allowed domain: $x")
                return nothing
            end
            return x
        end
    elseif mode == 1
        return (x, abs_sys) -> begin
            Xmap = SY.get_state_mapping(abs_sys)
            Xset = SY.get_state_domain(abs_sys)
            q = MP.get_state_by_coord(Xmap, x)
            if q !== nothing && MP.contains_state(Xset, Xmap, q)
                return x
            end
            # Otherwise: project to nearest allowed state (in coordinate space)
            # NOTE: this can be expensive if Xset is large.
            states = MP.enum_states(Xset, Xmap)

            isempty(states) && begin
                warn && @warn("Allowed state set is empty; cannot project $x")
                return nothing
            end

            # Choose a distance in full space or in selected dims
            function dist(qi::Int)
                xi = MP.get_coord_by_state(Xmap, qi)
                if dims === nothing
                    return LinearAlgebra.norm(xi - x)
                else
                    return LinearAlgebra.norm(xi[dims] - x[dims])
                end
            end

            qbest = argmin(dist, states)
            xproj = MP.get_coord_by_state(Xmap, qbest)

            warn && @warn(
                "State out of allowed domain: $x -> projected to state $qbest at $xproj (mode=1)"
            )
            return xproj
        end
    else
        error("Unknown mode=$mode")
    end
end

end # module
