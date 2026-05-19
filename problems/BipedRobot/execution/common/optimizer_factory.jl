include(joinpath(@__DIR__, "robot_setup.jl"))

using Dionysos
using MathematicalSystems
using StaticArrays
using LinearAlgebra
using JuMP
using MathOptInterface
using JLD2

const MOI = MathOptInterface
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

function save_optimizer(filename::AbstractString, optimizer)
    mkpath(dirname(filename))
    JLD2.jldopen(filename, "w") do file
        return file["optimizer"] = optimizer
    end
    return filename
end

function load_optimizer(filename::AbstractString)
    return JLD2.jldopen(filename, "r") do file
        return file["optimizer"]
    end
end

function build_robot_grids(concrete_system; simplify::Float64 = 3.0)
    n_state = MathematicalSystems.statedim(concrete_system)
    n_input = MathematicalSystems.inputdim(concrete_system)

    x0 = SVector{n_state, Float64}(zeros(n_state))
    hx =
        SVector{n_state, Float64}([fill(2π / 180, 3)..., fill(0.15, n_state - 3)...]) *
        simplify
    state_grid = MP.GridFree(x0, hx)

    u0 = SVector{n_input, Float64}(zeros(n_input))
    hu = SVector{n_input, Float64}(fill(1.0, n_input))
    input_grid = MP.GridFree(u0, hu)

    return state_grid, input_grid
end

function build_robot_problem(; tstep::Float64 = 0.1)
    return RobotProblem.problem(; robot_urdf = selected_robot_urdf(), tstep = tstep)
end

function build_robot_abstraction_optimizer(;
    execution_backend,
    simplify::Float64 = 3.0,
    tstep::Float64 = 0.1,
    state_filter = nothing,
    state_input_filter = nothing,
    print_level::Int = 2,
    progress_update_interval::Int = Int(1e2), # Update progress every x iterations (for printlevel >= 1)
)
    concrete_problem = build_robot_problem(; tstep = tstep)
    concrete_system = concrete_problem.system

    state_grid, input_grid = build_robot_grids(concrete_system; simplify = simplify)

    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_mapping"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("mapping_region"), concrete_system.X)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_filter"), state_filter)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_input_filter"), state_input_filter)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )

    MOI.set(optimizer, MOI.RawOptimizerAttribute("execution_backend"), execution_backend)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("progress_update_interval"),
        progress_update_interval,
    )

    return optimizer
end

function build_empty_abstraction_for_optimizer!(optimizer)
    abstraction_solver = optimizer.abstraction_solver
    abstraction_solver === nothing && error("Optimizer has no abstraction_solver.")

    return AB.UniformGridAbstraction.build_empty_abstraction!(abstraction_solver)
end
