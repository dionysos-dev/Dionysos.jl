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
import MathematicalSystems
MS = MathematicalSystems
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
    concrete_controller::Union{Nothing, MS.AbstractSystem, MS.AbstractMap}
    solve_time_sec::T
    print_level::Int
    randomize::Bool
    handle_out_of_domain::Function

    function Optimizer{T}() where {T}
        default_handler = (x, abs_sys) -> nothing  # stop if out of domain
        return new{T}(nothing, nothing, nothing, 0.0, 1, false, default_handler)
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

# Domain object whose membership is a predicate
struct PredicateDomain{F}
    pred::F
end
Base.in(x, X::PredicateDomain) = X.pred(x)

function solve_concrete_problem(
    abstract_system::Dionysos.Symbolic.GridBasedSymbolicModel,
    abstract_controller::MS.AbstractMap;
    randomize::Bool = false,
    handle_out_of_domain::Function = (x, abs_sys) -> nothing,  # default: stop
)
    k_abs = abstract_controller.h

    # map concrete x -> abstract state (or nothing)
    x_to_qs = function (x)
        state = SY.get_abstract_state(abstract_system, x)
        if state===nothing || !SY.is_allowed_state(abstract_system, state)
            xnew = handle_out_of_domain(x, abstract_system)
            xnew === nothing && return nothing
            statenew === nothing && return nothing
            SY.is_allowed_state(abstract_system, statenew) || return nothing
            return statenew
        end
        return state
    end

    # domain predicate on x (controller defined?)
    is_defined = function (x)
        qs = x_to_qs(x)
        qs === nothing && return false
        us = k_abs(qs)
        return us !== nothing && !(us isa AbstractVector && isempty(us))
    end

    # output map x -> concrete u (or nothing)
    f = function (x)
        qs = x_to_qs(x)
        qs === nothing && return nothing

        us = k_abs(qs)
        us === nothing && return nothing

        u_sym = if us isa AbstractVector
            isempty(us) ? nothing : (randomize ? rand(us) : first(us))
        else
            us
        end
        u_sym === nothing && return nothing

        return Dionysos.Symbolic.get_concrete_input(abstract_system, u_sym)
    end

    X = PredicateDomain(is_defined)
    nx::Int = Dionysos.Symbolic.get_state_dim(abstract_system)
    nu::Int = Dionysos.Symbolic.get_input_dim(abstract_system)
    return MS.ConstrainedBlackBoxMap(nx, nu, f, X)
end

# Concretizes an abstract finite-memory controller (on (qa,qs)) into a concrete
# finite-memory controller (on (qa,x)).
# Assumptions on `abstract_controller`:
# - output:  u_sym = h_abs(qa, qs)  (may be Int, Vector{Int}, or nothing)
# - update:  qa_next = g_abs(qa, qs_for_update)
# Returns:
# - concrete_controller :: MS.SystemWithOutput
# - qa0 is NOT stored here; you keep passing q0 to closed_loop_trajectory.
function solve_concrete_problem(
    abstract_system::Dionysos.Symbolic.GridBasedSymbolicModel,
    abstract_controller::MS.SystemWithOutput;
    randomize::Bool = false,
    handle_out_of_domain::Function = (x, abs_sys) -> nothing,
)
    # MS callables qa == abstract controller state
    h_abs = abstract_controller.outputmap.h   # (qa, qs) -> u_sym (or Vector) or nothing
    g_abs = MS.mapping(abstract_controller.s)         # (qa, qs) -> qa_next

    # concrete-state -> abstract-state (or nothing)
    x_to_qs = function (x)
        state = SY.get_abstract_state(abstract_system, x)
        if !SY.is_allowed_state(abstract_system, state)
            xnew = handle_out_of_domain(x, abstract_system)
            xnew === nothing && return nothing
            statenew = SY.get_abstract_state(abstract_system, xnew)
            !SY.is_allowed_state(abstract_system, statenew) || return nothing
            return statenew
        end
        return state
    end

    # output map for concrete controller: (qa, x) -> u_conc (or nothing)
    h_conc = function (qa, x)
        qs = x_to_qs(x)
        qs === nothing && return nothing

        us = h_abs((qa, qs))
        us === nothing && return nothing

        u_sym = if us isa AbstractVector
            isempty(us) ? nothing : (randomize ? rand(us) : first(us))
        else
            us
        end
        u_sym === nothing && return nothing

        return Dionysos.Symbolic.get_concrete_input(abstract_system, u_sym)
    end

    # state update for concrete controller: (qa, x_for_update) -> qa_next
    # (your rollout can choose x or x_next via update_on_next flag)
    g_conc = function (qa, x_for_update)
        qs = x_to_qs(x_for_update)
        qs === nothing && return qa     # or a dead-state if you prefer
        return g_abs(qa, qs)
    end

    # Domain predicate on (qa, x) for the outputmap
    is_defined_qax = function (qax)
        qa, x = qax
        qs = x_to_qs(x)
        qs === nothing && return false
        us = h_abs(qa, qs)
        return us !== nothing && !(us isa AbstractVector && isempty(us))
    end
    X_qax = PredicateDomain(is_defined_qax)

    # Build MS objects
    nx = Dionysos.Symbolic.get_state_dim(abstract_system)
    nu = Dionysos.Symbolic.get_input_dim(abstract_system)

    outmap = MS.ConstrainedBlackBoxMap(2, nu, qax -> begin
        qa, x = qax
        h_conc(qa, x)
    end, X_qax)

    memsys = MS.BlackBoxControlDiscreteSystem(
        (qa, x_for_update) -> g_conc(qa, x_for_update),
        1,
        nx,
    )

    return MS.SystemWithOutput(memsys, outmap)
end

function MOI.optimize!(optimizer::Optimizer)
    t_ref = time()

    # Ensure the concrete problem is defined
    if optimizer.abstraction_solver === nothing
        error("The concrete problem is not defined.")
    end

    # Compute abstraction if not already done
    if !is_abstraction_computed(optimizer)
        MOI.set(
            optimizer.abstraction_solver,
            MOI.RawOptimizerAttribute("print_level"),
            optimizer.print_level,
        )
        MOI.optimize!(optimizer.abstraction_solver)
    end

    # If there's a control solver, optimize it
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
        optimizer.concrete_controller = solve_concrete_problem(
            optimizer.abstraction_solver.abstract_system,
            abstract_controller;
            randomize = optimizer.randomize,
            handle_out_of_domain = optimizer.handle_out_of_domain,
        )
    end

    # Time elapsed
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
