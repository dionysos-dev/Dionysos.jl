export UniformGridAbstraction

module UniformGridAbstraction

import Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const SY = Dionysos.Symbolic
const DO = Dionysos.Domain

import StaticArrays: SVector, SMatrix
import MathematicalSystems
import HybridSystems
using JuMP

include("empty_problem.jl")
include("optimal_control_problem.jl")
include("safety_problem.jl")

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
    concrete_controller::Union{Nothing, ST.ContinuousController}
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

function solve_concrete_problem(
    abstract_system::Dionysos.Symbolic.GridBasedSymbolicModel,
    abstract_controller::Dionysos.System.SymbolicController,
)
    function is_defined(x)
        xpos = Dionysos.Domain.get_pos_by_coord(abstract_system.Xdom, x)
        if !(xpos ∈ abstract_system.Xdom)
            return false
        end
        source = Dionysos.Symbolic.get_state_by_xpos(abstract_system, xpos)
        if !Dionysos.System.is_defined(abstract_controller, source)
            return false
        end
        return true
    end
    function f(x; randomize::Bool = false)
        # 1. Abstract the concrete state
        xpos = Dionysos.Domain.get_pos_by_coord(abstract_system.Xdom, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("State out of domain: $x")
            return nothing
        end
        # 2. Get abstract state index
        source = Dionysos.Symbolic.get_state_by_xpos(abstract_system, xpos)
        # 3. Check if controller is defined there
        if !Dionysos.System.is_defined(abstract_controller, source)
            @warn "Uncontrollable state: $x"
            return nothing
        end
        # 4. Select a symbol from the admissible inputs
        inputs = Dionysos.System.get_all_controls(abstract_controller, source)
        symbol = randomize ? rand(inputs) : first(inputs)
        # 5. Map to concrete input
        u = Dionysos.Symbolic.get_concrete_input(abstract_system, symbol)

        return u
    end
    return ST.BlackBoxContinuousController(f, is_defined)
end

function is_abstraction_computed(optimizer::Optimizer)
    return optimizer.abstraction_solver !== nothing &&
           optimizer.abstraction_solver.abstract_system !== nothing
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

function get_input_symbol(controller, state)
    syms = Dionysos.Utils.fix_and_eliminate_first(controller, state)
    return isempty(syms) ? -1 : first(syms)[1]
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
