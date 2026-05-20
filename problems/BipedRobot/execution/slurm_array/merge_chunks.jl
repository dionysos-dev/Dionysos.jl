# ==============================================================================
# Merge SLURM-array transition chunks into one abstraction
# ==============================================================================

import Pkg
Pkg.instantiate()

using Printf

const SIMPLIFY = parse(Float64, get(ENV, "DIONYSOS_SIMPLIFY", "3.0"))
const TSTEP = parse(Float64, get(ENV, "DIONYSOS_TSTEP", "0.1"))

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

outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", default_abstraction_file())
nchunks = parse(Int, get(ENV, "DIONYSOS_NCHUNKS", "1"))

@info(
    "Preparing empty abstraction for merge",
    model = selected_robot_model(),
    nchunks = nchunks,
    outdir = outdir,
    outfile = outfile,
    simplify = SIMPLIFY,
    tstep = TSTEP,
)

execution_backend = SY.SlurmArrayBackend(nchunks, 1, outdir, :contiguous, true)

t_build_optimizer = @elapsed begin
    global optimizer = build_robot_abstraction_optimizer(;
        execution_backend = execution_backend,
        simplify = SIMPLIFY,
        tstep = TSTEP,
        print_level = 1,
    )
end

@printf("Optimizer build time:         %.3f s\n", t_build_optimizer)

abstract_system = build_empty_abstraction_for_optimizer!(optimizer)

@info "Merging transition chunks" nchunks outdir

t_merge = @elapsed SY.merge_transition_chunks!(
    abstract_system,
    outdir;
    nchunks = nchunks,
    print_level = 1,
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)

@printf("Merge time:                   %.3f s\n", t_merge)

t_save = @elapsed save_optimizer(outfile, optimizer)

@printf("Save time:                    %.3f s\n", t_save)

@info "Saved merged optimizer" outfile
