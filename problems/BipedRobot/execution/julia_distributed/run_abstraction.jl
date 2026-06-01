# ==============================================================================
# JuliaDistributed runner for robot abstraction
# ==============================================================================

import Pkg

USE_SYSIMAGE = false

DO_PKG_INSTANTIATE = !USE_SYSIMAGE
DO_PKG_PRECOMPILE = false

if DO_PKG_INSTANTIATE
    Pkg.instantiate()
end
if DO_PKG_PRECOMPILE
    Pkg.precompile()
end

using Distributed
using Printf

# ==============================================================================
# Parameters from environment
# ==============================================================================

const NWORKERS = parse(Int, get(ENV, "DIONYSOS_NWORKERS", "4"))

const USE_THREADED_PER_WORKER =
    lowercase(get(ENV, "DIONYSOS_THREADED_PER_WORKER", "false")) == "true"

const DISTRIBUTED_NPARTS = parse(Int, get(ENV, "DIONYSOS_NPARTS", string(NWORKERS)))

const PARTITION_STRATEGY = Symbol(get(ENV, "DIONYSOS_PARTITION_STRATEGY", "contiguous"))

tstep = 0.1

const BIPED_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const PROJECT_DIR = BIPED_ROOT

const SYSIMAGE_PATH = joinpath(PROJECT_DIR, "dionysos_robot_sysimage.dll")

# ==============================================================================
# Startup timing
# ==============================================================================

t_worker_warmup = 0.0

t_startup_total = @elapsed begin
    global t_master_packages = @elapsed begin
        using Dionysos
        using MathematicalSystems
        using StaticArrays
        using LinearAlgebra
        using JuMP
        using MathOptInterface
        using JLD2
    end

    global const MOI = MathOptInterface
    global const DI = Dionysos
    global const UT = DI.Utils
    global const ST = DI.System
    global const MP = DI.Mapping
    global const SY = DI.Symbolic
    global const OP = DI.Optim
    global const AB = OP.Abstraction

    global t_master_includes = @elapsed begin
        include(joinpath(@__DIR__, "..", "common", "robot_setup.jl"))

        global robot_problem_path = selected_robot_problem_path()
        include(robot_problem_path)
        using .RobotProblem
    end

    global t_worker_creation = @elapsed begin
        n_to_add = max(NWORKERS - nworkers(), 0)

        if n_to_add > 0
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

    global t_worker_warmup = @elapsed begin
        if !isempty(workers())
            @info "Warming up workers..." nworkers = length(workers())

            futures = [
                remotecall(
                    RobotProblem.warmup_robot_problem!,
                    p;
                    robot_urdf = selected_robot_urdf(),
                    tstep = tstep,
                ) for p in workers()
            ]

            fetch.(futures)

            @info "Worker warm-up finished."
        end
    end
end

@printf("Startup total time:           %.3f s\n", t_startup_total)
@printf("  Master package load:        %.3f s\n", t_master_packages)
@printf("  Master file includes:       %.3f s\n", t_master_includes)
@printf("  Worker creation:            %.3f s\n", t_worker_creation)
@printf("  Worker package load:        %.3f s\n", t_worker_packages)
@printf("  Worker file includes:       %.3f s\n", t_worker_includes)
@printf("  Worker warm-up:             %.3f s\n", t_worker_warmup)

println("Workers available: ", length(workers()))
println("USE_SYSIMAGE: ", USE_SYSIMAGE)

# ==============================================================================
# Optimizer factory
# ==============================================================================

include(joinpath(@__DIR__, "..", "common", "optimizer_factory.jl"))

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

robot_urdf = selected_robot_urdf()
tstep = 0.1
domain = RobotProblem.default_robot_domain()

simplify = 3.0
discretization = RobotDiscretizationConfig(;
    hx = SVector(2π / 180, 2π / 180, 2π / 180, 0.15, 0.15, 0.15) * simplify,
    hu = SVector(1.0, 1.0, 1.0) * simplify,
)

# ==============================================================================
# Abstraction
# ==============================================================================

@info(
    "Starting Julia distributed abstraction",
    nworkers = length(workers()),
    nparts = DISTRIBUTED_NPARTS,
    partition_strategy = PARTITION_STRATEGY,
    threaded_per_worker = USE_THREADED_PER_WORKER,
)

execution_backend = SY.JuliaDistributedBackend(
    nothing,                  # use Distributed.workers()
    DISTRIBUTED_NPARTS,
    PARTITION_STRATEGY,
    USE_THREADED_PER_WORKER,
    true,                     # warmup abstraction workers
)

concrete_problem = RobotProblem.problem(;
    robot_urdf = robot_urdf,
    tstep = tstep,
    domain = domain,
    Δt_simu = 1e-4,
    simulator = :history,
)

optimizer = build_robot_abstraction_optimizer(
    concrete_problem,
    execution_backend,
    discretization;
    print_level = 2,
    progress_update_interval = Int(1e3),
    save_concrete_traj = false,
)

t_opt_wall = @elapsed MOI.optimize!(optimizer)

t_construct =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))

@printf("Abstraction optimize! wall time: %.3f s\n", t_opt_wall)
@printf("Abstraction reported construction time: %.3f s\n", t_construct)

outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", default_abstraction_file())

t_save = @elapsed save_optimizer(outfile, optimizer)

@printf("Abstraction save time: %.3f s\n", t_save)
@info "Saved abstraction" outfile
