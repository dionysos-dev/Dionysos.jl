# ==============================================================================
#  Runner script for abstraction + optimal-control simulations
# ==============================================================================

using MathematicalSystems
using StaticArrays
using LinearAlgebra
using Plots
using JuMP
using JLD2

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(@__DIR__, "..", "src", "RS_tools.jl"))
import .RS_tools

include(joinpath(@__DIR__, "robot_problem.jl"))
include(joinpath(@__DIR__, "utils.jl"))

# ==============================================================================
# Script parameters
# ==============================================================================
const FILENAME = joinpath(@__DIR__, "Abstraction.jld2")

const COMPUTE_ABSTRACTION = false
const SAVE_ABSTRACTION = true
const LOAD_ABSTRACTION = false

const SIMULATE_FIRST_STEP = true
const SIMULATE_SECOND_STEP = false

robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf")
tstep = 0.1

# ==============================================================================
# Helpers
# ==============================================================================
reached_target(problem) = (x -> (x ∈ problem.target_set))

function build_optimizer(; concrete_problem, state_grid, input_grid)
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_update_interval"), Int(1e3))
    MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_dt"), 60)

    return optimizer
end

function solve_and_simulate!(
    optimizer,
    concrete_system,
    xstart::SVector,
    target_low::SVector,
    target_high::SVector;
    nstep::Int = 300,
    out_of_domain_handler = nothing,
)
    # Problem
    I = UT.HyperRectangle(xstart, xstart)   # start forced in cell of xstart
    T = UT.HyperRectangle(target_low, target_high)

    problem = DI.Problem.OptimalControlProblem(
        concrete_system,
        I,
        T,
        nothing,
        nothing,
        DI.Problem.Infinity(),
    )

    # Solve abstract problem
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    if out_of_domain_handler !== nothing
        MOI.set(
            optimizer,
            MOI.RawOptimizerAttribute("handle_out_of_domain"),
            out_of_domain_handler,
        )
    end

    MOI.optimize!(optimizer)

    t_abs = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
    println("Time to solve the abstract problem: $(t_abs) sec")

    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # Simulate closed-loop
    stopfun = reached_target(problem)
    x_traj, u_traj = ST.get_closed_loop_trajectory(
        concrete_system,
        controller,
        xstart,
        nstep;
        stopping = stopfun,
    )

    return x_traj, u_traj
end

function make_test_trajectory(; N = 300, dt = 0.1)
    seq = Vector{SVector{6, Float64}}(undef, N)

    for k in 1:N
        t = (k - 1) * dt

        # Joint angles (rad)
        q3 = 5π/180 * sin(0.5t)
        q4 = 4π/180 * sin(0.5t + π/4)
        q5 = 6π/180 * sin(0.5t + π/2)

        # Joint velocities (rad/s)
        v3 = 5π/180 * 0.5 * cos(0.5t)
        v4 = 4π/180 * 0.5 * cos(0.5t + π/4)
        v5 = 6π/180 * 0.5 * cos(0.5t + π/2)

        seq[k] = @SVector [q3, q4, q5, v3, v4, v5]
    end

    return ST.Trajectory(seq)
end

# ==============================================================================
# System setup
# ==============================================================================
concrete_problem = RobotProblem.problem(; robot_urdf = robot_urdf, tstep = tstep)
concrete_system = concrete_problem.system

n_state = MathematicalSystems.statedim(concrete_system)
n_input = MathematicalSystems.inputdim(concrete_system)
println("n_state: ", n_state)
println("n_input: ", n_input)

# ==============================================================================
# Abstraction (compute / save / load)
# ==============================================================================
optimizer = nothing

if COMPUTE_ABSTRACTION
    x0 = SVector{n_state, Float64}(zeros(n_state))
    hx = SVector{n_state, Float64}([fill(2π/180, 3)..., fill(0.15, 3)...])
    state_grid = DO.GridFree(x0, hx)

    u0 = SVector{n_input, Float64}(zeros(n_input))
    hu = SVector{n_input, Float64}(fill(1.0, n_input))
    input_grid = DO.GridFree(u0, hu)

    optimizer = build_optimizer(;
        concrete_problem = concrete_problem,
        state_grid = state_grid,
        input_grid = input_grid,
    )

    MOI.optimize!(optimizer)

    if SAVE_ABSTRACTION
        AB.UniformGridAbstraction.export_abstraction_jld2(optimizer, FILENAME)
        println("Saved abstraction to: ", FILENAME)
    end
end

if LOAD_ABSTRACTION
    optimizer = AB.UniformGridAbstraction.load_abstraction_jld2(FILENAME)
    println("Loaded abstraction from: ", FILENAME)
end

# ==============================================================================
# Simulations
# ==============================================================================
if SIMULATE_FIRST_STEP
    println("\nFirst step:\n")
    x0 = SVector{n_state, Float64}(zeros(n_state))

    t_low = SVector{n_state, Float64}([-12π/180, 7π/180, 8π/180, -0.75, -0.30, -0.30])
    t_high = SVector{n_state, Float64}([-8π/180, 9π/180, 12π/180, 0.30, 0.75, 0.75])

    x_traj, u_traj =
        solve_and_simulate!(optimizer, concrete_system, x0, t_low, t_high; nstep = 300)
    # x_traj = make_test_trajectory()

    println(x_traj, "\n")
    println(u_traj, "\n")
    rs, vis = RS_tools.get_visualization_tool(; robot_urdf = robot_urdf)
    RS_tools.animate_trajectory!(vis, x_traj.seq; dt = tstep)
end

if SIMULATE_SECOND_STEP
    println("\nSecond step:\n")

    x0 = SVector{n_state, Float64}([
        -0.15352800685754736,
        0.11944498327439435,
        0.21311298746900986,
        0.0,
        0.0,
        0.0,
    ])

    t_low =
        SVector{n_state, Float64}([-1.1π/180, -1.1π/180, -1.1π/180, -0.75, -0.30, -0.30])
    t_high = SVector{n_state, Float64}([1.1π/180, 1.1π/180, 1.1π/180, 0.30, 0.75, 0.75])

    handler = AB.UniformGridAbstraction.make_out_of_domain_handler(; mode = 1, warn = true)

    x_traj, u_traj = solve_and_simulate!(
        optimizer,
        concrete_system,
        x0,
        t_low,
        t_high;
        nstep = 300,
        out_of_domain_handler = handler,
    )

    println(x_traj, "\n")
    println(u_traj, "\n")
end
