# ==============================================================================
# Merge SLURM-array transition chunks into one abstraction
# ==============================================================================

import Pkg
Pkg.instantiate()

using Printf

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

simplify = 1.0 # 3.0
discretization = RobotDiscretizationConfig(;
    hx = SVector(2π / 180, 2π / 180, 2π / 180, 0.15, 0.15, 0.15) * simplify,
    hu = SVector(1.0, 1.0, 1.0) * simplify,
)

outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", default_abstraction_file())
nchunks = 200

# ------------------------------------------------------------------------------
# Merge
# ------------------------------------------------------------------------------

@info(
    "Preparing empty abstraction for merge",
    nchunks = nchunks,
    outdir = outdir,
    outfile = outfile,
)

execution_backend = SY.SlurmArrayBackend(nchunks, 1, outdir, :contiguous, true)

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

abstract_system = build_empty_abstraction_for_optimizer!(optimizer)

@info "Merging transition chunks" nchunks outdir

t_merge = @elapsed SY.merge_transition_chunks!(
    abstract_system,
    outdir;
    nchunks = nchunks,
    print_level = 1,
)

println(SY.has_metadata(abstract_system))
tr = first(SY.enum_transitions(abstract_system))
println(tr)
println(SY.get_metadata(abstract_system, tr))

MOI.set(optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)

@printf("Merge time:                   %.3f s\n", t_merge)

t_save = @elapsed save_optimizer(outfile, optimizer)

@printf("Save time:                    %.3f s\n", t_save)

@info "Saved merged optimizer" outfile
