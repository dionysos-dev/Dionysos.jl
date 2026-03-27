# ==============================================================================
#  Runner script for abstraction + optimal-control simulations
# ==============================================================================

# ---------------------------------------------------------------------------
#  Configuration from environment variables
# ---------------------------------------------------------------------------
const _t_script_start = time()

const USE_DISTRIBUTED = lowercase(get(ENV, "DIONYSOS_DISTRIBUTED", "false")) == "true"
const USE_THREADED = lowercase(get(ENV, "DIONYSOS_THREADED", "false")) == "true"
const N_PARTS = parse(Int, get(ENV, "DIONYSOS_NPARTS", "300"))
const N_PROCS = N_PARTS

using Distributed
if USE_DISTRIBUTED && length(workers()) < 2
    addprocs(max(N_PROCS, 2) - length(workers()))
end

using MathematicalSystems
using StaticArrays
using LinearAlgebra
using JuMP
using JLD2

@everywhere using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

robot_problem_path = joinpath(@__DIR__, "robot_problem.jl")
@everywhere include($robot_problem_path)

@everywhere include(joinpath(@__DIR__, "utils.jl"))

const loading_time = time() - _t_script_start
println("TIMING loading_time: ", loading_time)

# ==============================================================================
# Script parameters
# ==============================================================================
const MODE_TAG =
    USE_DISTRIBUTED && USE_THREADED ? "hybrid" :
    USE_DISTRIBUTED ? "distributed" : USE_THREADED ? "threaded" : "serial"
const OUTDIR = get(ENV, "DIONYSOS_OUTDIR", @__DIR__)
mkpath(OUTDIR)
const FILENAME = joinpath(OUTDIR, "Abstraction_$(MODE_TAG).jld2")

const COMPUTE_ABSTRACTION = true
const SAVE_ABSTRACTION = false
const LOAD_ABSTRACTION = false

const SIMULATE_FIRST_STEP = false
const SIMULATE_SECOND_STEP = false

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

    MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed"), USE_DISTRIBUTED)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed_nparts"), N_PARTS)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("distributed_partition_strategy"),
        :roundrobin,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), USE_THREADED)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)

    # MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_update_interval"), Int(1e5))
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
    I = UT.HyperRectangle(xstart, xstart)
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
robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf")
tstep = 0.1

const _t_system_setup_start = time()
concrete_problem = RobotProblem.problem(; robot_urdf = robot_urdf, tstep = tstep)
concrete_system = concrete_problem.system
const system_setup_time = time() - _t_system_setup_start
println("TIMING system_setup_time: ", system_setup_time)

n_state = MathematicalSystems.statedim(concrete_system)
n_input = MathematicalSystems.inputdim(concrete_system)

println("n_state: ", n_state)
println("n_input: ", n_input)

const pre_abstraction_time = time() - _t_script_start
println("TIMING pre_abstraction_time: ", pre_abstraction_time)

state_filter = nothing # RobotProblem.in_gait_tube
state_input_filter = nothing # RobotProblem.input_allowed

# ==============================================================================
# Abstraction (compute / save / load)
# ==============================================================================
optimizer = nothing

if COMPUTE_ABSTRACTION
    x0 = SVector{n_state, Float64}(zeros(n_state))
    hx = SVector{n_state, Float64}([fill(2π/180, 3)..., fill(0.15, 3)...])*2.5
    state_grid = MP.GridFree(x0, hx)

    u0 = SVector{n_input, Float64}(zeros(n_input))
    hu = SVector{n_input, Float64}(fill(1.0, n_input))
    input_grid = MP.GridFree(u0, hu)

    optimizer = build_optimizer(;
        concrete_problem = concrete_problem,
        state_grid = state_grid,
        input_grid = input_grid,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )

    RobotProblem.reset_rbd_timing!()
    MOI.optimize!(optimizer)

    abstraction_time =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
    println("Time to construct the abstraction: $(abstraction_time)")

    # Collect RBD timing: on master for serial/threaded, from workers for distributed
    if USE_DISTRIBUTED && length(workers()) > 1
        local rbd_total_ns = Int64(0)
        local rbd_total_calls = Int64(0)
        for wid in workers()
            wstats = @fetchfrom wid RobotProblem.get_rbd_timing()
            rbd_total_ns += Int64(round(wstats.time_sec * 1e9))
            rbd_total_calls += wstats.call_count
            println(
                "TIMING rbd_worker worker_id=$wid time=$(wstats.time_sec) calls=$(wstats.call_count)",
            )
        end
        println("TIMING rbd_simulate_time: ", rbd_total_ns / 1e9)
        println("TIMING rbd_call_count: ", rbd_total_calls)
    else
        rbd_stats = RobotProblem.get_rbd_timing()
        println("TIMING rbd_simulate_time: ", rbd_stats.time_sec)
        println("TIMING rbd_call_count: ", rbd_stats.call_count)
    end

    global post_abstraction_time = time() - _t_script_start
    println("TIMING post_abstraction_time: ", post_abstraction_time)

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
if SIMULATE_FIRST_STEP
    println("\nFirst step:\n")
    x0 = SVector{n_state, Float64}(zeros(n_state))

    t_low = SVector{n_state, Float64}([-12π/180, 7π/180, 8π/180, -0.75, -0.30, -0.30])
    t_high = SVector{n_state, Float64}([-8π/180, 9π/180, 12π/180, 0.30, 0.75, 0.75])

    x_traj, u_traj =
        solve_and_simulate!(optimizer, concrete_system, x0, t_low, t_high; nstep = 300)

    # Save trajectory to file (works on headless servers)
    traj_path = joinpath(OUTDIR, "trajectory_step1_$(MODE_TAG).jld2")
    @save traj_path x_traj u_traj
    println("Trajectory saved to: ", traj_path)
    println("Final state: ", x_traj.seq[end])
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

    # Save trajectory to file (works on headless servers)
    traj_path = joinpath(OUTDIR, "trajectory_step2_$(MODE_TAG).jld2")
    @save traj_path x_traj u_traj
    println("Trajectory saved to: ", traj_path)
    println("Final state: ", x_traj.seq[end])
end

rmprocs(workers())
println("Workers removed")
