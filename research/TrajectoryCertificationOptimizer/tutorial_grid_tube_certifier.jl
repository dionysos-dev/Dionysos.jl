# Tutorial 1 — the uniform-grid TUBE certifier on the 2-D integrator.
#
# Same modular pipeline as tutorial 0, but the certification block is the
# uniform-grid trajectory certifier: a grid tube is inflated around the
# candidate trajectory, a local abstraction is built inside it, and the
# certificate is the grid controllable set driving every tube cell into the
# target — set-valued (cells), where the ellipsoidal certifier is funnel-valued
# (ellipsoids + feedback gains).
#
#     julia --project=test research/TrajectoryCertificationOptimizer/tutorial_grid_tube_certifier.jl
#
# Companions: tutorial_pipeline_blocks.jl (the pipeline block by block),
# tutorial_ellipsoidal_knobs.jl (the ellipsoidal LMI knobs, annotated).

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

import LazySets
import LinearAlgebra as LA
import MathOptInterface as MOI
import MathematicalSystems as MS
using StaticArrays
using Random
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
state_grid = MP.GridFree(SVector(-2.0, -2.0), SVector(0.2, 0.2))
input_grid = MP.GridFree(SVector(-1.0, -1.0), SVector(0.5, 0.5))

# ------------------------------------------------------------
# 1) Composite generator: abstraction seed → MPPI polish
# ------------------------------------------------------------

abstraction_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    PR.AlternatingSimulationProblem(concrete_system, concrete_system.X),
)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
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

mppi_generator = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = nothing,             # the composite block wires the seed in
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

combo_gen = AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
    abstraction_generator,
    mppi_generator;
    Δt = Δt,
    num_substeps = 5,
)

# ------------------------------------------------------------
# 2) Uniform-grid tube certifier
# ------------------------------------------------------------

const UGTC = AB.UniformGridTrajectoryCertifier

local_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(local_optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(local_optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(local_optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    local_optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(
    local_optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    Integrator.jacobian_bound(),
)
MOI.set(local_optimizer, MOI.Silent(), true)

ug_certifier = UGTC.TrajectoryCertifier(;
    optimizer = local_optimizer,
    radius = 0.6,                  # tube half-width around the candidate
    margin = 0.0,
    incl_mode = MP.INNER,
    n_between = 0,
    max_step = nothing,
    enforce_safe_max_step = true,
    handle_system_domain = true,
)

# ------------------------------------------------------------
# 3) The generate ⇄ certify driver
# ------------------------------------------------------------

println("— trajectory-certification driver (grid tube) —")
tc_optimizer = AB.TrajectoryCertificationOptimizer.Optimizer(combo_gen, ug_certifier)
MOI.set(tc_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.optimize!(tc_optimizer)

composite_traj = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("trajectory"))
controller = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("controller"))
println(
    "  success = $(MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("success"))), ",
    "status = $(MOI.get(tc_optimizer, MOI.TerminationStatus())), ",
    "time = $(round(MOI.get(tc_optimizer, MOI.SolveTimeSec()); digits = 2)) s, ",
    "controller = $(typeof(controller))",
)

# The certificate lives on the grid: the controllable set of the local
# abstraction built inside the tube.
ug_optimizer = tc_optimizer.trajectory_certifier.optimizer
abstract_system = MOI.get(ug_optimizer, MOI.RawOptimizerAttribute("abstract_system"))
controllable_set = MOI.get(ug_optimizer, MOI.RawOptimizerAttribute("controllable_set"))
XMapping = SY.get_state_mapping(abstract_system)

# ------------------------------------------------------------
# 4) Static plot + dashboard animation
# ------------------------------------------------------------

fig = plot(; aspect_ratio = :equal, title = "Grid-tube certificate on the 2-D integrator")
plot!(fig, concrete_problem; aspect_ratio = :equal)

if composite_traj !== nothing
    X_tube = UGTC.build_tube(
        composite_traj,
        ug_certifier.radius;
        margin = ug_certifier.margin,
        n_between = ug_certifier.n_between,
        max_step = ug_certifier.max_step,
        enforce_safe_max_step = ug_certifier.enforce_safe_max_step,
        X_domain = ug_certifier.handle_system_domain ? concrete_problem.system.X : nothing,
    )
    plot!(fig, X_tube; alpha = 0.25, color = :blue, label = "certified tube")
    plot!(
        (controllable_set, XMapping);
        color = :yellow,
        linecolor = :yellow,
        label = "controllable set",
    )
    plot!(
        fig,
        composite_traj;
        color = :green,
        dims = [1, 2],
        linewidth = 2,
        label = "certified trajectory",
    )
end
plot_path = joinpath(@__DIR__, "tutorial_grid_tube_certifier.png")
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
    composite_traj;
    xdims = (1, 2),
    udims = (1, 2),
    Δt = Δt,
    fps = 10,
    filename = joinpath(@__DIR__, "tutorial_grid_tube_certifier_dashboard.gif"),
    title = "2-D integrator, grid-tube certified",
)
println("— dashboard saved: $dash_path —")
