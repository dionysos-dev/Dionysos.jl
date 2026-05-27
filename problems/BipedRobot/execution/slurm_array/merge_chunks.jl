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
# # Cedric & Maxime settings
# domain = RobotProblem.default_robot_domain()
# simplify = 1.0 # 3.0
# discretization = RobotDiscretizationConfig(;
#     hx = SVector(2π / 180, 2π / 180, 2π / 180, 0.15, 0.15, 0.15) * simplify,
#     hu = SVector(1.0, 1.0, 1.0) * simplify,
# )
# nchunks = 200

# Baptiste settings
domain = RobotProblem.RobotDomainConfig(;
    x_lb = SVector{6, Float64}([-25*π/180, -25*π/180, -10*π/180, -2, -1, -2.5]),
    x_ub = SVector{6, Float64}([25*π/180, 25*π/180, 80*π/180, 1, 2, 2.5]),
    u_lb = SVector{3, Float64}((-4.0, -4.0, -4.0)),
    u_ub = SVector{3, Float64}((4.0, 4.0, 4.0)),
)
discretization = RobotDiscretizationConfig(;
    hx = SVector{6, Float64}([
        0.025490685402655037,
        0.021796907174709057,
        0.025842360657025935,
        0.6553157054436949,
        0.4732745368241606,
        1.548368319358735,
    ]),
    hu = SVector{3, Float64}((1.3221879606782365, 1.7200722827025006, 0.7845631194998841)),
)
nchunks = 5000

outdir = get(ENV, "DIONYSOS_TRANSITION_OUTDIR", default_transition_outdir())
outfile = get(ENV, "DIONYSOS_ABSTRACTION_FILE", default_abstraction_file())

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
