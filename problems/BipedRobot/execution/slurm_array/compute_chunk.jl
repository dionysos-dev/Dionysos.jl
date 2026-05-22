# ==============================================================================
# SLURM-array chunk computation for robot abstraction
# ==============================================================================

import Pkg
Pkg.instantiate()

using Printf

# ------------------------------------------------------------------------------
# Startup timing
# ------------------------------------------------------------------------------

t_startup_total = @elapsed begin
    global t_packages = @elapsed begin
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

    global t_includes = @elapsed begin
        include(joinpath(@__DIR__, "..", "common", "robot_setup.jl"))

        global robot_problem_path = selected_robot_problem_path()
        include(robot_problem_path)
        using .RobotProblem
    end
end

@printf("Startup total time:           %.3f s\n", t_startup_total)
@printf("  Package load:               %.3f s\n", t_packages)
@printf("  File includes:              %.3f s\n", t_includes)

include(joinpath(@__DIR__, "..", "common", "optimizer_factory.jl"))

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

robot_urdf = selected_robot_urdf()
tstep = 0.1
domain = RobotProblem.default_robot_domain()
# RobotDomainConfig(
#     x_lb = SVector(-0.3, 0.0, 0.0, -0.4, -0.2, -0.2),
#     x_ub = SVector(0.0, 0.3, 0.4, 0.2, 0.4, 0.4),
#     u_lb = SVector(-2.0, -2.0, -3.0),
#     u_ub = SVector(2.0, 2.0, 3.0),
# )

simplify = 1.0 # 3.0
discretization = RobotDiscretizationConfig(;
    hx = SVector(2π / 180, 2π / 180, 2π / 180, 0.15, 0.15, 0.15) * simplify,
    hu = SVector(1.0, 1.0, 1.0) * simplify,
)

# ------------------------------------------------------------------------------
# Chunk parameters
# ------------------------------------------------------------------------------

outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
nchunks = parse(Int, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))
chunk_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

@info(
    "Starting SLURM abstraction chunk",
    chunk_id = chunk_id,
    nchunks = nchunks,
    outdir = outdir,
)

execution_backend = SY.SlurmArrayBackend(
    nchunks,
    chunk_id,
    outdir,
    :contiguous,
    true, # write_only: do not add transitions to local abstract system
)

global t_warmup = @elapsed begin
    @info "Warming up robot problem"
    RobotProblem.warmup_robot_problem!(;
        robot_urdf = robot_urdf,
        tstep = tstep,
        Δt_simu = 5e-4,
        simulator = :custom,
    )
    @info "Robot warm-up finished"
end
@printf("  Robot warm-up:              %.3f s\n", t_warmup)

# ------------------------------------------------------------------------------
# Build and compute chunk
# ------------------------------------------------------------------------------

concrete_problem = RobotProblem.problem(;
    robot_urdf = robot_urdf,
    tstep = tstep,
    domain = domain,
    Δt_simu = 5e-4, # 1e-4,
    simulator = :custom,
)

t_build_optimizer = @elapsed begin
    global optimizer = build_robot_abstraction_optimizer(
        concrete_problem,
        execution_backend,
        discretization;
        print_level = 2,
        progress_update_interval = Int(1e3),
        save_concrete_traj = true,
    )
end

@printf("Optimizer build time:         %.3f s\n", t_build_optimizer)

t_opt_wall = @elapsed MOI.optimize!(optimizer)

@printf("Chunk optimize! wall time:    %.3f s\n", t_opt_wall)

@info(
    "Finished SLURM abstraction chunk",
    chunk_id = chunk_id,
    nchunks = nchunks,
    elapsed = t_opt_wall,
)
