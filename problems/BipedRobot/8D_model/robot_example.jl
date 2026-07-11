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
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(@__DIR__, "..", "src", "RS_tools.jl"))
import .RS_tools

include(joinpath(@__DIR__, "robot_problem.jl"))

# ==============================================================================
# Script parameters
# ==============================================================================
const FILENAME = joinpath(@__DIR__, "Abstraction.jld2")

const COMPUTE_ABSTRACTION = true
const SAVE_ABSTRACTION = false
const LOAD_ABSTRACTION = false

const SIMULATE = false

# ==============================================================================
# Helpers
# ==============================================================================
reached_target(problem) = (x -> (x ∈ problem.target_set))

function build_optimizer(;
    concrete_problem,
    state_grid,
    input_grid,
    state_filter = nothing,
    state_input_filter = nothing,
)
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_filter"), state_filter)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_input_filter"), state_input_filter)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )

    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.Silent(), true)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("execution_backend"),
        SY.SequentialBackend(),
    )
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
    I = UT.box(xstart, xstart)   # start forced in cell of xstart
    T = UT.box(target_low, target_high)

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

function make_test_trajectory_8d(; N::Int = 300, dt::Real = 0.1)
    seq = Vector{SVector{8, Float64}}(undef, N)

    # amplitudes (rad) and a common frequency (rad/s)
    A = @SVector [5π/180, 4π/180, 6π/180, 3π/180]  # LH, RH, LK, RK
    ω = 0.5
    ϕ = @SVector [0.0, π/4, π/2, 3π/4]

    for k in 1:N
        t = (k - 1) * dt

        # positions
        LH = A[1] * sin(ω*t + ϕ[1])
        RH = A[2] * sin(ω*t + ϕ[2])
        LK = A[3] * sin(ω*t + ϕ[3])
        RK = A[4] * sin(ω*t + ϕ[4])

        # velocities (derivative)
        dLH = A[1] * ω * cos(ω*t + ϕ[1])
        dRH = A[2] * ω * cos(ω*t + ϕ[2])
        dLK = A[3] * ω * cos(ω*t + ϕ[3])
        dRK = A[4] * ω * cos(ω*t + ϕ[4])

        seq[k] = @SVector [LH, RH, LK, RK, dLH, dRH, dLK, dRK]
    end

    return ST.Trajectory(seq)
end

# ==============================================================================
# System setup
# ==============================================================================

robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf")
tstep = 0.1

concrete_problem = RobotProblem.problem(; robot_urdf = robot_urdf, tstep = tstep)
concrete_system = concrete_problem.system

n_state = MathematicalSystems.statedim(concrete_system)
n_input = MathematicalSystems.inputdim(concrete_system)
println("n_state: ", n_state)
println("n_input: ", n_input)

state_filter = nothing # RobotProblem.in_gait_tube
state_input_filter = nothing # RobotProblem.input_allowed

# ==============================================================================
# Abstraction (compute / save / load)
# ==============================================================================
optimizer = nothing

if COMPUTE_ABSTRACTION
    x0 = SVector{n_state, Float64}(zeros(n_state))
    hx = SVector{n_state, Float64}(fill(0.3, n_state))
    state_grid = MP.GridFree(x0, hx)

    u0 = SVector{n_input, Float64}(zeros(n_input))
    hu = SVector{n_input, Float64}(fill(3.0, n_input))
    input_grid = MP.GridFree(u0, hu)

    optimizer = build_optimizer(;
        concrete_problem = concrete_problem,
        state_grid = state_grid,
        input_grid = input_grid,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )

    MOI.optimize!(optimizer)

    if SAVE_ABSTRACTION
        AB.UniformGridAbstraction.export_abstraction_jld2(optimizer, FILENAME)
        println("Saved abstraction to: ", FILENAME)
    end
end

if LOAD_ABSTRACTION
    optimizer = AB.UniformGridAbstraction.import_abstraction_jld2(FILENAME)
    println("Loaded abstraction from: ", FILENAME)
end

# ==============================================================================
# Simulations
# ==============================================================================
if SIMULATE
    println("\nFirst step:\n")
    x0 = SVector{n_state, Float64}(zeros(n_state))

    t_low = SVector{n_state, Float64}([0.1, 0.1, -0.1, -0.1, -0.8, -0.8, -0.8, -0.8])
    t_high = SVector{n_state, Float64}([0.5, 0.5, -0.5, -0.5, 0.8, 0.8, 0.8, 0.8])

    x_traj, u_traj =
        solve_and_simulate!(optimizer, concrete_system, x0, t_low, t_high; nstep = 300)
    # x_traj = make_test_trajectory()

    rs, vis = RS_tools.get_visualization_tool(; robot_urdf = robot_urdf)
    RS_tools.animate_trajectory!(vis, x_traj.seq; dt = tstep)
end
