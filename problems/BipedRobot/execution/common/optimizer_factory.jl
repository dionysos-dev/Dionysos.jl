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

Base.@kwdef struct RobotDiscretizationConfig{X, U}
    hx::X
    hu::U
end

function default_robot_discretization(; scale::Float64 = 1.0)
    return RobotDiscretizationConfig(;
        hx = scale .* SVector(2π / 180, 2π / 180, 2π / 180, 0.15, 0.15, 0.15),
        hu = scale .* SVector(1.0, 1.0, 1.0),
    )
end

function build_robot_grids(concrete_system, discretization)
    n_state = MathematicalSystems.statedim(concrete_system)
    n_input = MathematicalSystems.inputdim(concrete_system)

    length(discretization.hx) == n_state ||
        error("Expected hx of length $n_state, got $(length(discretization.hx)).")

    length(discretization.hu) == n_input ||
        error("Expected hu of length $n_input, got $(length(discretization.hu)).")

    x0 = SVector{n_state, Float64}(zeros(n_state))
    u0 = SVector{n_input, Float64}(zeros(n_input))

    hx = SVector{n_state, Float64}(discretization.hx)
    hu = SVector{n_input, Float64}(discretization.hu)

    state_grid = MP.GridFree(x0, hx)
    input_grid = MP.GridFree(u0, hu)

    return state_grid, input_grid
end

function build_robot_abstraction_optimizer(
    concrete_problem,
    execution_backend,
    discretization;
    state_filter = nothing,
    state_input_filter = nothing,
    print_level::Int = 2,
    progress_update_interval::Int = Int(1e2),
    save_concrete_traj = false,
)
    concrete_system = concrete_problem.system

    state_grid, input_grid = build_robot_grids(concrete_system, discretization)

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
    if save_concrete_traj
        MOI.set(
            optimizer,
            MOI.RawOptimizerAttribute("transition_metadata"),
            SY.TransitionMetadata(),
        ) # SY.NoTransitionMetadata(), SY.TransitionMetadata()
    end

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
