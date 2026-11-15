export UniformGridAbstraction

module UniformGridAbstraction

import Dionysos
const UT = Dionysos.Utils
const PR = Dionysos.Problem
const ST = Dionysos.System
const SY = Dionysos.Symbolic
const DO = Dionysos.Domain

import StaticArrays: SVector, SMatrix
import MathematicalSystems
MS = MathematicalSystems
import HybridSystems
import Spot
using JuMP

export Optimizer

include("empty_problem.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")
include("cosafe_LTL_problem.jl")

"""
    Optimizer{T} <: MOI.AbstractOptimizer

A high-level abstraction-based solver that automatically orchestrates system abstraction and control synthesis.  
This wrapper follows the **classical abstraction pipeline** (e.g., as in SCOTS), where the state and input spaces are discretized into hyper-rectangular cells, independent of the specific control task.

It delegates responsibility to modular sub-solvers: one for abstraction and one for control, depending on the type of problem to be solved.

---

### Structure and Sub-solvers

The optimizer internally manages two sub-solvers:

- `abstraction_solver`: [`OptimizerEmptyProblem`](@ref):  
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
    abstraction_solver::Union{Nothing, OptimizerEmptyProblem{T}}
    control_solver::Union{Nothing, MOI.AbstractOptimizer}
    concrete_controller::Union{Nothing, MS.AbstractSystem, MS.AbstractMap}
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
        model.abstraction_solver = OptimizerEmptyProblem()
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
        if isa(value, Dionysos.Problem.EmptyProblem)
            model.abstraction_solver = OptimizerEmptyProblem()
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("empty_problem"),
                value,
            )
            model.control_solver = nothing  # Pas de solveur de contrôle
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
        if model.abstraction_solver.empty_problem === nothing
            empty_problem = Dionysos.Problem.EmptyProblem(value.system, value.system.X)
            MOI.set(
                model.abstraction_solver,
                MOI.RawOptimizerAttribute("empty_problem"),
                empty_problem,
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

# Domain object whose membership is a predicate
struct PredicateDomain{F}
    pred::F
end
Base.in(x, X::PredicateDomain) = X.pred(x)

function solve_concrete_problem(
    abstract_system::Dionysos.Symbolic.GridBasedSymbolicModel,
    abstract_controller::MS.AbstractMap;
    randomize::Bool = false,
)
    k_abs = abstract_controller.h

    # concrete-state -> abstract-state (or nothing)
    x_to_qs = function (x)
        xpos = Dionysos.Domain.get_pos_by_coord(abstract_system.Xdom, x)
        (xpos ∈ abstract_system.Xdom) || return nothing
        return Dionysos.Symbolic.get_state_by_xpos(abstract_system, xpos)
    end

    # implicit domain predicate on x
    is_defined = function (x)
        qs = x_to_qs(x)
        qs === nothing && return false
        us = k_abs(qs)
        return us !== nothing && !(us isa AbstractVector && isempty(us))
    end

    # output map x -> u (or nothing)
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
    nx::Int = Dionysos.Symbolic.get_concrete_state_dim(abstract_system)
    nu::Int = Dionysos.Symbolic.get_concrete_input_dim(abstract_system)
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
)
    # MS callables qa == abstract controller state
    h_abs = abstract_controller.outputmap.h   # (qa, qs) -> u_sym (or Vector) or nothing
    g_abs = MS.mapping(abstract_controller.s)         # (qa, qs) -> qa_next

    # concrete-state -> abstract-state (or nothing)
    x_to_qs = function (x)
        xpos = Dionysos.Domain.get_pos_by_coord(abstract_system.Xdom, x)
        (xpos ∈ abstract_system.Xdom) || return nothing
        return Dionysos.Symbolic.get_state_by_xpos(abstract_system, xpos)
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
    nx = Dionysos.Symbolic.get_concrete_state_dim(abstract_system)
    nu = Dionysos.Symbolic.get_concrete_input_dim(abstract_system)

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
            abstract_controller,
        )
    end

    # Time elapsed
    optimizer.solve_time_sec = time() - t_ref
    return
end

using DataFrames, CSV

function export_controller_csv(
    optimizer::UniformGridAbstraction.Optimizer,
    filename::String,
)
    abstract_system = optimizer.abstraction_solver.abstract_system
    abstract_controller = optimizer.control_solver.abstract_controller
    abstract_controller === nothing && error("Controller not available")

    return export_controller_csv(abstract_system, abstract_controller, filename)
end

function export_controller_csv(abstract_system, abstract_controller, basename::String)
    grid = SY.get_state_grid(abstract_system)

    CSV.write(basename * "_Grid.csv", build_grid_df(grid); delim = ';')
    CSV.write(
        basename * "_StateMap.csv",
        build_state_map_df(abstract_system, grid);
        delim = ';',
    )
    CSV.write(
        basename * "_ControllerMap.csv",
        build_controller_map_df(abstract_system, abstract_controller);
        delim = ';',
    )
    return CSV.write(
        basename * "_InputMap.csv",
        build_input_map_df(abstract_system);
        delim = ';',
    )
end

function build_grid_df(grid)
    origin = DO.get_origin(grid)
    h = DO.get_h(grid)
    ndims = length(origin)

    header = ["key"; ["x$(j)" for j in 1:ndims]]
    rows = [["origin"; origin], ["h"; h]]

    df = DataFrame()
    for j in 1:length(header)
        df[!, Symbol(header[j])] = getindex.(rows, j)
    end
    return df
end

function build_state_map_df(abstract_system, grid)
    ndims = length(DO.get_h(grid))
    headers = ["abstract_state"; ["x$(j)" for j in 1:ndims]]
    states = SY.enum_states(abstract_system)
    rows = [(s, SY.get_xpos_by_state(abstract_system, s)...) for s in states]
    return DataFrame([headers[i] => getindex.(rows, i) for i in 1:length(headers)])
end

function build_controller_map_df(abstract_system, abstract_controller)
    states = SY.enum_states(abstract_system)
    rows = [(s, get_input_symbol(abstract_controller, s)) for s in states]
    return DataFrame([
        "abstract_state" => getindex.(rows, 1),
        "abstract_input" => getindex.(rows, 2),
    ])
end

function get_input_symbol(controller, state; randomize = false)
    !(state in controller.X) && return -1

    u = controller.h(state)

    u === nothing && return -1
    (u isa AbstractVector && isempty(u)) && return -1

    return u isa AbstractVector ? (randomize ? rand(u) : first(u)) : u
end

function build_input_map_df(abstract_system)
    inputs = SY.enum_inputs(abstract_system)
    ndims_u = length(SY.get_concrete_input(abstract_system, first(inputs)))
    headers = ["abstract_input"; ["u$(j)" for j in 1:ndims_u]]
    rows = [(i, SY.get_concrete_input(abstract_system, i)...) for i in inputs]
    return DataFrame([headers[i] => getindex.(rows, i) for i in 1:length(headers)])
end

function load_controller_data_csv(basename::String)
    grid_df = CSV.read(basename * "_Grid.csv", DataFrame; delim = ';')
    state_df = CSV.read(basename * "_StateMap.csv", DataFrame; delim = ';')
    ctrl_df = CSV.read(basename * "_ControllerMap.csv", DataFrame; delim = ';')
    input_df = CSV.read(basename * "_InputMap.csv", DataFrame; delim = ';')
    return parse_controller_tables(grid_df, state_df, ctrl_df, input_df)
end

function parse_controller_tables(grid_df, state_df, ctrl_df, input_df)
    origin = Vector{Float64}(grid_df[grid_df.key .== "origin", Not(:key)][1, :])
    h = Vector{Float64}(grid_df[grid_df.key .== "h", Not(:key)][1, :])

    pos2state = Dict{Vector{Int}, Int}()
    for row in eachrow(state_df)
        pos = [Float64(row[Symbol("x$i")]) for i in 1:(ncol(state_df) - 1)]
        pos2state[pos] = row.abstract_state
    end

    state2input = Dict(ctrl_df.abstract_state .=> ctrl_df.abstract_input)

    input2u = Dict{Int, Vector{Float64}}()
    for row in eachrow(input_df)
        u = [Float64(row[Symbol("u$i")]) for i in 1:(ncol(input_df) - 1)]
        input2u[row.abstract_input] = u
    end

    return origin, h, pos2state, state2input, input2u
end

end
