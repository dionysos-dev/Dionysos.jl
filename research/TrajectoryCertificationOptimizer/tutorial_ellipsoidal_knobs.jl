# Tutorial 2 — the ellipsoidal certifier's knobs, annotated, on the 2-D
# integrator.
#
# The toy system is exactly affine: the local affine approximation is exact,
# the Jacobians are constant, the Lipschitz remainder is identically zero, and
# the nonlinear error bounds vanish. Consequently `linearization_δx/δu` do not
# affect the model at all — their only role is through the box-consistency
# checks of the certification algorithm. Any failure here can be attributed to
# LMI feasibility, numerical conditioning, input constraints, or the
# consistency checks — never to linearization error. That makes this file the
# right place to see what each ChainOptions knob does.
#
#     julia --project=test research/TrajectoryCertificationOptimizer/tutorial_ellipsoidal_knobs.jl
#
# Companions: tutorial_pipeline_blocks.jl (the pipeline block by block),
# tutorial_grid_tube_certifier.jl (uniform-grid TUBE certifier).

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

_X_ = LazySets.Hyperrectangle(; low = SVector(-3.0, -3.0), high = SVector(4.0, 4.0))
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
    AB.UniformGridAbstraction.CENTER_SIMULATION,   # cheap: the dynamics are exact
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
# 2) Ellipsoidal backward certifier — every knob spelled out
# ------------------------------------------------------------

Symbolics.@variables x[1:2] u[1:2] w[1:2]

# Symbolic discrete map the certifier linearizes (exactly: x⁺ = x + Δt·(u + w)).
fsymbolic = [x[1] + Δt * (u[1] + w[1]), x[2] + Δt * (u[2] + w[2])]

# Zero disturbance set — the LMI API still expects a noise format.
Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))

# The provider tells the certifier how to build local affine approximations
# (and Lipschitz remainder bounds) along the trajectory.
provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,                    # symbolic dynamics
    collect(x),                   # symbolic state variables
    collect(u),                   # symbolic input variables
    collect(w),                   # symbolic disturbance variables
    [0.0, 0.0],                   # disturbance radius ΔW
    ST.format_input_set(_U_),     # input constraints for the LMI
    ST.format_noise_set(Wset),    # disturbance vertices for the LMI
)

# Optional quadratic transition cost: the LMIs then also bound the worst-case
# cost-to-go of the synthesized feedback (here: pure input energy).
trans_cost = UT.QuadraticStateControlFunction(
    zeros(2, 2),                  # Q: state cost
    Matrix{Float64}(LA.I, 2, 2),  # R: input cost
    zeros(2, 2),                  # N
    zeros(2),                     # q
    zeros(2),                     # r
    0.0,                          # v
)

ellip_opts = EB.ChainOptions(;
    # Upper bounds on the predecessor-ellipsoid size and on the controller's
    # input deviation from the nominal (RADII at this boundary; the kernels
    # bound the squared deviations). Larger makes the LMI easier.
    maxδx = 1e6,
    maxδu = 1e6,
    # Objective tradeoff: min λ·cost − (1−λ)·size. Small λ favors LARGER
    # ellipsoids; large λ favors lower transition cost.
    λ = 0.3,
    # Shape matrix Q of the terminal ellipsoid, centered at the trajectory
    # endpoint. `nothing` would inscribe it in the target box instead.
    terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.5^2,
    transition_cost = trans_cost,
    # Fixed linearization boxes (exact dynamics ⟹ only the consistency checks
    # see them: a certified ellipsoid/controller image larger than its box is
    # rejected by the box-consistency gate, not by any approximation error).
    linearization_δx = [1.1, 1.0],
    linearization_δu = [0.5, 0.5],
    # Size objective: :maximin maximizes the smallest semi-axis (collapse-proof
    # default); :logdet maximizes true volume (fragile on hard chains — C2
    # measured 2/8); :trace is the linear volume proxy.
    objective = :maximin,
)

sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)
certifier = EB.BackwardCertifier(provider, sdp, ellip_opts)

# ------------------------------------------------------------
# 3) The generate ⇄ certify driver
# ------------------------------------------------------------

println("— trajectory-certification driver (ellipsoidal) —")
tc_optimizer = AB.TrajectoryCertificationOptimizer.Optimizer(combo_gen, certifier)
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

cert_result = EB.get_result(certifier)
println(
    "  chain: $(length(cert_result.lmi_data.ellipsoids)) ellipsoids, ",
    "success = $(cert_result.success), failed_k = $(cert_result.failed_k)",
)

# Per-step records carry the LMI diagnostics — the first place to look when a
# chain dies (status, adaptive-box outcome, gate reasons).
for s in cert_result.steps
    s.status == :ok || println("  k=$(s.k) status=$(s.status)")
end

# ------------------------------------------------------------
# 4) Static plot: state funnels + the controller's input-ellipsoid images
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

# κ(Eₖ) — the image of a funnel ellipsoid under its affine feedback. The
# controller is NOT forced to interpolate the nominal input: if the nominal is
# near the boundary of U, the optimizer shifts the certified controller inward
# so the whole image stays feasible.
function input_image_ellipsoid(E, κ)
    K = Matrix{Float64}(κ.A)
    b = collect(Float64, κ.c)
    c_u = K * collect(Float64, LazySets.center(E)) + b
    Q_u = K * Matrix{Float64}(LazySets.shape_matrix(E)) * K'
    return LazySets.Ellipsoid(c_u, UT._symmetrize(Q_u))
end

fig = plot(; aspect_ratio = :equal, title = "Ellipsoidal certificate, annotated")
plot!(fig, concrete_problem; aspect_ratio = :equal)
plot_ellipsoid_chain!(fig, cert_result; max_ellipsoids = 100)
if composite_traj !== nothing
    plot!(
        fig,
        composite_traj;
        color = :green,
        dims = [1, 2],
        label = "certified trajectory",
    )
end

fig_u = plot(; aspect_ratio = :equal, title = "Input set and κ(Eₖ) images")
plot!(fig_u, _U_; alpha = 0.2, label = "input set U")
valid_steps = filter(s -> s.ellipsoid !== nothing && s.kappa !== nothing, cert_result.steps)
for (j, s) in enumerate(valid_steps[1:max(1, div(length(valid_steps), 20)):end])
    plot!(
        fig_u,
        input_image_ellipsoid(s.ellipsoid, s.kappa);
        linewidth = 1.5,
        linestyle = :dash,
        alpha = 0.6,
        label = j == 1 ? "κ(Eₖ)" : "",
    )
end

final_fig = plot(fig, fig_u; layout = (1, 2), size = (1400, 620))
plot_path = joinpath(@__DIR__, "tutorial_ellipsoidal_knobs.png")
savefig(final_fig, plot_path)
println("— plot saved: $plot_path —")

# ------------------------------------------------------------
# 5) Dashboard animation
# ------------------------------------------------------------

sys_plot = Integrator.system_plot!(; xlims = (-3.2, 4.2), ylims = (-3.2, 4.2))
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
    filename = joinpath(@__DIR__, "tutorial_ellipsoidal_knobs_dashboard.gif"),
    title = "2-D integrator, ellipsoidal certificate",
)
println("— dashboard saved: $dash_path —")
