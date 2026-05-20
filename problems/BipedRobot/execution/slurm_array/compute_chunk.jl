# ==============================================================================
# SLURM-array chunk computation for robot abstraction
# ==============================================================================

import Pkg
Pkg.instantiate()

using Printf

const SIMPLIFY = parse(Float64, get(ENV, "DIONYSOS_SIMPLIFY", "3.0"))
const TSTEP = parse(Float64, get(ENV, "DIONYSOS_TSTEP", "0.1"))

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

    global t_warmup = @elapsed begin
        if isdefined(RobotProblem, :warmup_robot_problem!)
            @info "Warming up robot problem" robot_urdf = selected_robot_urdf() tstep =
                TSTEP
            RobotProblem.warmup_robot_problem!(;
                robot_urdf = selected_robot_urdf(),
                tstep = TSTEP,
            )
            @info "Robot warm-up finished"
        end
    end
end

@printf("Startup total time:           %.3f s\n", t_startup_total)
@printf("  Package load:               %.3f s\n", t_packages)
@printf("  File includes:              %.3f s\n", t_includes)
@printf("  Robot warm-up:              %.3f s\n", t_warmup)

include(joinpath(@__DIR__, "..", "common", "optimizer_factory.jl"))

# ------------------------------------------------------------------------------
# Chunk parameters
# ------------------------------------------------------------------------------

outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
nchunks = parse(Int, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))
chunk_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

@info(
    "Starting SLURM abstraction chunk",
    model = selected_robot_model(),
    chunk_id = chunk_id,
    nchunks = nchunks,
    outdir = outdir,
    simplify = SIMPLIFY,
    tstep = TSTEP,
)

execution_backend = SY.SlurmArrayBackend(
    nchunks,
    chunk_id,
    outdir,
    :contiguous,
    true,      # write_only: do not add transitions to local abstract system
)

# ------------------------------------------------------------------------------
# Build and compute chunk
# ------------------------------------------------------------------------------

t_build_optimizer = @elapsed begin
    global optimizer = build_robot_abstraction_optimizer(;
        execution_backend = execution_backend,
        simplify = SIMPLIFY,
        tstep = TSTEP,
        print_level = 2,
        progress_update_interval = Int(1e2),
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
