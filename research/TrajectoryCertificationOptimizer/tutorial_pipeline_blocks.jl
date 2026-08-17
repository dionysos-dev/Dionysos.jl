# Tutorial 0 — the modular pipeline, block by block, on the 2-D integrator.
#
# Every block of the trajectory-certification pipeline is instantiated and run
# on the simplest possible system (ẋ = u, exactly affine — zero linearization
# error), so each interface can be inspected in isolation:
#
#   1) abstraction-based trajectory generator (uniform grid + Dijkstra),
#   2) MPPI refinement (typed stage costs, fixed rng),
#   3) composite generator (abstraction seed → MPPI polish),
#   4) ellipsoidal backward certifier, standalone,
#   5) the generate ⇄ certify driver tying them together,
#   6) static plot + dashboard animation.
#
#     julia --project=test research/TrajectoryCertificationOptimizer/tutorial_pipeline_blocks.jl
#
# Companions: tutorial_grid_tube_certifier.jl (uniform-grid TUBE certifier),
# tutorial_ellipsoidal_knobs.jl (the ellipsoidal LMI knobs, annotated).

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const AB = DI.Optim.Abstraction
const EB = AB.EllipsoidalTrajectoryCertifier

import LazySets
import LinearAlgebra as LA
import MathOptInterface as MOI
import MathematicalSystems as MS
using StaticArrays
using Random
using Symbolics
import Clarabel
using JuMP: optimizer_with_attributes
using Plots

const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")
include(joinpath(PROBLEMS, "Integrator", "integrator.jl"))

# ------------------------------------------------------------
# 0) Problem: box domain, two-component target, no obstacles
# ------------------------------------------------------------

_X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(4.0, 4.0))
_U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
_I_ = LazySets.Hyperrectangle(; low = SVector(-1.7, -1.7), high = SVector(-1.6, -1.6))
g11 = LazySets.Hyperrectangle(; low = SVector(-1.0, 3.0), high = SVector(-0.3, 3.7))
g12 = LazySets.Hyperrectangle(; low = SVector(1.0, 2.0), high = SVector(3.0, 3.7))
target_set = UT.set_union([g11, g12])

concrete_system = Integrator.system(; _X_ = _X_, _U_ = _U_)
concrete_problem =
    PR.OptimalControlProblem(concrete_system, _I_, target_set, nothing, nothing, 20)

Δt = 0.3

# ------------------------------------------------------------
# 1) Trajectory generator block 1: abstraction-based (uniform grid)
# ------------------------------------------------------------

println("— abstraction generator —")
abstraction_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    PR.AlternatingSimulationProblem(concrete_system, concrete_system.X),
)
MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("state_grid"),
    MP.GridFree(SVector(-2.0, -2.0), SVector(0.2, 0.2)),
)
MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("input_grid"),
    MP.GridFree(SVector(-1.0, -1.0), SVector(0.5, 0.5)),
)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    Integrator.jacobian_bound(),
)
MOI.set(abstraction_optimizer, MOI.Silent(), true)

abstraction_generator = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
    abstraction_optimizer;
    initial_state = SVector(-1.65, -1.65),
    concrete = false,
    nstep = 50,
)
AB.set_problem!(abstraction_generator, concrete_problem)
AB.generate!(abstraction_generator)
println(
    "  success = $(AB.get_success(abstraction_generator)), ",
    "time = $(round(AB.get_solve_time(abstraction_generator); digits = 2)) s",
)
abstraction_traj = AB.get_trajectory(abstraction_generator)

# ------------------------------------------------------------
# 2) Trajectory generator block 2: MPPI refinement (typed stage costs)
# ------------------------------------------------------------

println("— MPPI —")
discrete_problem = PR.discretize_problem(concrete_problem, Δt)

mppi_generator = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = abstraction_traj,
    nstep = 50,
    nsamples = 1000,
    niter = 20,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.3, 0.3)),
    project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0)),
    cost = AB.CompositeCost(
        AB.ReachObjectiveCost(target_set),
        AB.InputEffortCost(0.01),
        AB.DomainPenaltyCost(concrete_system.X),
    ),
)
AB.set_problem!(mppi_generator, discrete_problem)
AB.generate!(mppi_generator)
println(
    "  success = $(AB.get_success(mppi_generator)), ",
    "time = $(round(AB.get_solve_time(mppi_generator); digits = 2)) s",
)
mppi_traj = AB.get_trajectory(mppi_generator)

# ------------------------------------------------------------
# 3) Composite generator: abstraction seed → MPPI polish, one block
# ------------------------------------------------------------

println("— composite generator —")
combo_gen = AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
    abstraction_generator,
    mppi_generator;
    Δt = Δt,
    num_substeps = 5,
)
AB.set_problem!(combo_gen, concrete_problem)
AB.generate!(combo_gen)
composite_traj = AB.get_trajectory(combo_gen)
composite_seed = AB.get_seed(combo_gen)
println("  success = $(AB.get_success(combo_gen))")

# ------------------------------------------------------------
# 4) Ellipsoidal backward certifier, standalone
# ------------------------------------------------------------

println("— ellipsoidal certifier (standalone) —")
Symbolics.@variables x[1:2] u[1:2] w[1:2]
fsymbolic = [x[1] + Δt * (u[1] + w[1]), x[2] + Δt * (u[2] + w[2])]
Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))
provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    collect(x),
    collect(u),
    collect(w),
    [0.0, 0.0],
    ST.format_input_set(_U_),
    ST.format_noise_set(Wset),
)

# Fixed linearization boxes: the integrator is exactly affine, so the boxes
# carry no approximation content — only the box-consistency checks see them.
ellip_opts = EB.ChainOptions(;
    maxδx = 30.0,
    maxδu = 1.0,
    λ = 0.05,
    terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.5^2,
    linearization_δx = [0.2, 0.2],
    linearization_δu = [0.5, 0.5],
    objective = :maximin,          # :logdet (volume) and :trace are the alternatives
)

sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)
certifier = EB.BackwardCertifier(provider, sdp, ellip_opts)
AB.set_problem!(certifier, concrete_problem)
AB.set_trajectory!(certifier, composite_traj)
AB.certify!(certifier)
println(
    "  success = $(AB.get_success(certifier)), ",
    "time = $(round(AB.get_solve_time(certifier); digits = 2)) s",
)
cert_result = EB.get_result(certifier)

# ------------------------------------------------------------
# 5) The generate ⇄ certify driver: same blocks, one optimizer
# ------------------------------------------------------------

println("— trajectory-certification driver —")
tc_optimizer = AB.TrajectoryCertificationOptimizer.Optimizer(
    AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
        abstraction_generator,
        mppi_generator;
        Δt = Δt,
        num_substeps = 5,
    ),
    EB.BackwardCertifier(provider, sdp, ellip_opts),
)
MOI.set(tc_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.optimize!(tc_optimizer)
driver_traj = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("trajectory"))
controller = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("controller"))
println(
    "  success = $(MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("success"))), ",
    "status = $(MOI.get(tc_optimizer, MOI.TerminationStatus())), ",
    "time = $(round(MOI.get(tc_optimizer, MOI.SolveTimeSec()); digits = 2)) s, ",
    "controller = $(typeof(controller))",
)

# ------------------------------------------------------------
# 6) Static plot + dashboard animation
# ------------------------------------------------------------

function plot_ellipsoid_chain!(fig, cert_result; max_ellipsoids = 100)
    cert_result === nothing && return fig
    ellipsoids = cert_result.lmi_data.ellipsoids
    isempty(ellipsoids) && return fig
    idxs = unique(
        round.(
            Int,
            range(1, length(ellipsoids); length = min(max_ellipsoids, length(ellipsoids))),
        ),
    )
    for (j, idx) in enumerate(idxs)
        plot!(
            fig,
            ellipsoids[idx];
            label = j == 1 ? "backward funnel" : "",
            linewidth = 1.5,
            linestyle = :dash,
            alpha = 0.7,
        )
    end
    return fig
end

fig = plot(; aspect_ratio = :equal, title = "Modular pipeline on the 2-D integrator")
plot!(fig, concrete_problem; aspect_ratio = :equal)
plot_ellipsoid_chain!(fig, cert_result; max_ellipsoids = 10)
plot!(fig, abstraction_traj; color = :blue, dims = [1, 2], label = "abstraction trajectory")
plot!(fig, mppi_traj; color = :red, dims = [1, 2], label = "MPPI trajectory")
plot!(fig, composite_seed; color = :orange, dims = [1, 2], label = "composite seed")
plot!(fig, driver_traj; color = :green, dims = [1, 2], label = "certified trajectory")
plot_path = joinpath(@__DIR__, "tutorial_pipeline_blocks.png")
savefig(fig, plot_path)
println("— plot saved: $plot_path —")

sys_plot = Integrator.system_plot!(; xlims = (-2.2, 4.2), ylims = (-2.2, 4.2))
dash_plot! = (f, xk, uk) -> begin
    plot!(f, _I_; color = :gray, alpha = 0.4)
    plot!(f, g11; color = :green, alpha = 0.4)
    plot!(f, g12; color = :green, alpha = 0.4)
    sys_plot(f, xk, uk)
end
dash_path = DI.animate_trajectory_dashboard(
    dash_plot!,
    driver_traj;
    xdims = (1, 2),
    udims = (1, 2),
    Δt = Δt,
    fps = 10,
    filename = joinpath(@__DIR__, "tutorial_pipeline_blocks_dashboard.gif"),
    title = "2-D integrator, certified",
)
println("— dashboard saved: $dash_path —")
