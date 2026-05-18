# ==============================================================================
# Runner script for abstraction + optimal-control simulations
# ==============================================================================

import Pkg

USE_SYSIMAGE = false
NWORKERS = 2 # parse(Int, get(ENV, "DIONYSOS_NWORKERS", "1"))

# Only do package management in non-sysimage / setup mode if desired
DO_PKG_INSTANTIATE = !USE_SYSIMAGE
DO_PKG_PRECOMPILE = false
PROJECT_DIR = abspath(joinpath(@__DIR__, ".."))

if DO_PKG_INSTANTIATE
    Pkg.activate(PROJECT_DIR)
    Pkg.instantiate()
end
if DO_PKG_PRECOMPILE
    Pkg.precompile()
end

using Distributed
using Printf

robot_problem_path = joinpath(@__DIR__, "robot_problem.jl")
utils_path = joinpath(@__DIR__, "utils.jl")
rsviz_path = joinpath(@__DIR__, "..", "src", "RSVisualization.jl")

# project used by workers
PROJECT_DIR = abspath(joinpath(@__DIR__, ".."))

# sysimage built in problems/BipedRobot
SYSIMAGE_PATH = joinpath(PROJECT_DIR, "dionysos_robot_sysimage.dll")

# ------------------------------------------------------------------------------
# Timed startup
# ------------------------------------------------------------------------------
t_startup_total = @elapsed begin
    global t_master_packages = @elapsed begin
        using Dionysos
        using MathematicalSystems
        using StaticArrays
        using LinearAlgebra
        using JuMP
        using Plots
        using JLD2
        using MathOptInterface
    end

    global const MOI = MathOptInterface
    global const DI = Dionysos
    global const UT = DI.Utils
    global const ST = DI.System
    global const MP = DI.Mapping
    global const OP = DI.Optim
    global const AB = OP.Abstraction

    global t_master_includes = @elapsed begin
        include(robot_problem_path)
        include(utils_path)
        using .RobotProblem
    end

    global t_worker_creation = @elapsed begin
        if length(workers()) < NWORKERS
            n_to_add = NWORKERS - length(workers())

            if USE_SYSIMAGE
                addprocs(
                    n_to_add;
                    exeflags = `--project=$(PROJECT_DIR) --sysimage=$(SYSIMAGE_PATH)`,
                )
            else
                addprocs(n_to_add; exeflags = `--project=$(PROJECT_DIR)`)
            end
        end
    end

    global t_worker_packages = @elapsed begin
        @everywhere begin
            using Dionysos
        end
    end

    global t_worker_includes = @elapsed begin
        @everywhere include($robot_problem_path)
        @everywhere using .RobotProblem
    end
end

@printf("Startup total time:           %.3f s\n", t_startup_total)
@printf("  Master package load:        %.3f s\n", t_master_packages)
@printf("  Master file includes:       %.3f s\n", t_master_includes)
@printf("  Worker creation:            %.3f s\n", t_worker_creation)
@printf("  Worker package load:        %.3f s\n", t_worker_packages)
@printf("  Worker file includes:       %.3f s\n", t_worker_includes)
println("Workers available: ", length(workers()))
println("USE_SYSIMAGE: ", USE_SYSIMAGE)

# ==============================================================================
# Script parameters
# ==============================================================================
FILENAME = joinpath(@__DIR__, "Abstraction.jld2")

COMPUTE_ABSTRACTION = true
SAVE_ABSTRACTION = false
LOAD_ABSTRACTION = false

SIMULATE_FIRST_STEP = false
SIMULATE_SECOND_STEP = false

USE_DISTRIBUTED = length(workers()) > 0
USE_THREADED_PER_WORKER = false
DISTRIBUTED_NPARTS = length(workers())
DISTRIBUTED_PARTITION_STRATEGY = :contiguous # :roundrobin, :contiguous
SIMPLIFY = 1.5 # increase to simplify abstraction (e.g. by increasing grid size)
println("Simplify : ", SIMPLIFY)

# Only load visualization tools on master, and only if needed
if SIMULATE_FIRST_STEP || SIMULATE_SECOND_STEP
    include(rsviz_path)
    using .RSVisualization
end

# ==============================================================================
# Helpers
# ==============================================================================
reached_target(problem) = (x -> (x ∈ problem.target_set))

function warmup_workers!(; robot_urdf, tstep)
    isempty(workers()) && return nothing
    @info "Warming up workers..." nworkers = length(workers())

    futures = [
        remotecall(
            RobotProblem.warmup_robot_problem!,
            p;
            robot_urdf = robot_urdf,
            tstep = tstep,
        ) for p in workers()
    ]

    fetch.(futures)

    @info "Worker warm-up finished."
    return nothing
end

function build_optimizer(;
    concrete_problem,
    state_grid,
    input_grid,
    state_filter = nothing,
    state_input_filter = nothing,
    distributed = false,
    distributed_nparts = 1,
    distributed_partition_strategy = :roundrobin,
    threaded = false,
    print_level = 2,
)
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    concrete_system = concrete_problem.system

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

    MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed"), distributed)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed_nparts"), distributed_nparts)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("distributed_partition_strategy"),
        distributed_partition_strategy,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), threaded)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_update_interval"), Int(1e2))
    MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_dt"), 60.0)

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

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)

    if out_of_domain_handler !== nothing
        MOI.set(
            optimizer,
            MOI.RawOptimizerAttribute("handle_out_of_domain"),
            out_of_domain_handler,
        )
    end

    t_solve_wall = @elapsed MOI.optimize!(optimizer)
    @printf("Abstract problem wall time:     %.3f s\n", t_solve_wall)

    controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

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

        q3 = 5π / 180 * sin(0.5 * t)
        q4 = 4π / 180 * sin(0.5 * t + π / 4)
        q5 = 6π / 180 * sin(0.5 * t + π / 2)

        v3 = 5π / 180 * 0.5 * cos(0.5 * t)
        v4 = 4π / 180 * 0.5 * cos(0.5 * t + π / 4)
        v5 = 6π / 180 * 0.5 * cos(0.5 * t + π / 2)

        seq[k] = @SVector [q3, q4, q5, v3, v4, v5]
    end

    return ST.Trajectory(seq)
end

# ==============================================================================
# System setup
# ==============================================================================
robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf")
tstep = 0.1

t_global_setup = time()

t_warmup = 0.0
if USE_DISTRIBUTED
    t_warmup = @elapsed warmup_workers!(; robot_urdf = robot_urdf, tstep = tstep)
    @printf("Worker warm-up time: %.3f s\n", t_warmup)
end

t_problem_build = @elapsed begin
    global concrete_problem = RobotProblem.problem(; robot_urdf = robot_urdf, tstep = tstep)
    global concrete_system = concrete_problem.system
end
@printf("Concrete problem construction time: %.3f s\n", t_problem_build)

n_state = MathematicalSystems.statedim(concrete_system)
n_input = MathematicalSystems.inputdim(concrete_system)

println("n_state: ", n_state)
println("n_input: ", n_input)
println("distributed: ", USE_DISTRIBUTED)
println("nworkers: ", length(workers()))
println("distributed_nparts: ", DISTRIBUTED_NPARTS)
println("threaded_per_worker: ", USE_THREADED_PER_WORKER)

state_filter = nothing
state_input_filter = nothing
# state_filter = RobotProblem.in_gait_tube
# state_input_filter = RobotProblem.input_allowed

# ==============================================================================
# Abstraction
# ==============================================================================
optimizer = nothing

if COMPUTE_ABSTRACTION
    x0 = SVector{n_state, Float64}(zeros(n_state))
    hx = SVector{n_state, Float64}([fill(2π / 180, 3)..., fill(0.15, 3)...]) * SIMPLIFY
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
        distributed = USE_DISTRIBUTED,
        distributed_nparts = DISTRIBUTED_NPARTS,
        distributed_partition_strategy = DISTRIBUTED_PARTITION_STRATEGY,
        threaded = USE_THREADED_PER_WORKER,
        print_level = 2,
    )

    t_opt_wall = @elapsed MOI.optimize!(optimizer)
    @printf("Abstraction optimize! wall time: %.3f s\n", t_opt_wall)

    t_construct =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
    @printf("Abstraction reported construction time: %.3f s\n", t_construct)

    if SAVE_ABSTRACTION
        t_save =
            @elapsed AB.UniformGridAbstraction.export_abstraction_jld2(optimizer, FILENAME)
        @printf("Abstraction save time: %.3f s\n", t_save)
        println("Saved abstraction to: ", FILENAME)
    end
end

if LOAD_ABSTRACTION
    t_load = @elapsed global optimizer =
        AB.UniformGridAbstraction.import_abstraction_jld2(FILENAME)
    @printf("Abstraction load time: %.3f s\n", t_load)
    println("Loaded abstraction from: ", FILENAME)
end

@printf("Total setup-to-abstraction time so far: %.3f s\n", time() - t_global_setup)

# ==============================================================================
# Simulations
# ==============================================================================
if SIMULATE_FIRST_STEP
    println("\nFirst step:\n")
    x0 = SVector{n_state, Float64}(zeros(n_state))

    t_low = SVector{n_state, Float64}([-12π / 180, 7π / 180, 8π / 180, -0.75, -0.30, -0.30])
    t_high = SVector{n_state, Float64}([-8π / 180, 9π / 180, 12π / 180, 0.30, 0.75, 0.75])

    # x_traj, u_traj =
    #     solve_and_simulate!(optimizer, concrete_system, x0, t_low, t_high; nstep = 300)
    x_traj = make_test_trajectory(; N = 300, dt = tstep)

    rs, vis = RSVisualization.get_visualization_tool(; robot_urdf = robot_urdf)
    RSVisualization.animate_trajectory!(vis, x_traj.seq; dt = tstep)
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

    t_low = SVector{n_state, Float64}([
        -1.1π / 180,
        -1.1π / 180,
        -1.1π / 180,
        -0.75,
        -0.30,
        -0.30,
    ])
    t_high =
        SVector{n_state, Float64}([1.1π / 180, 1.1π / 180, 1.1π / 180, 0.30, 0.75, 0.75])

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
