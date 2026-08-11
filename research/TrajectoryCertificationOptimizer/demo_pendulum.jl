# Certified pendulum swing-up — the P6 pendulum demo (plan.md §8.1).
#
# The modern pipeline end-to-end: abstraction seed → MPPI (stage costs,
# ESS-adaptive temperature, antithetic sampling, input reserve) → periodic unwrap
# → ellipsoidal certification (maximin objective, adaptive boxes, default
# inscribed terminal) in both directions + the bidirectional handoff. The
# convex objective ("reachability_up_medium_power_no_obstacle") has a plain box
# input set, so the certifier's input set IS the plant's input set — no ±10.5
# fiction (plan.md §4.2-E).
#
#     julia --project=bench research/TrajectoryCertificationOptimizer/demo_pendulum.jl
#
# Status (2026-08-11, Clarabel, fixed rng): FULLY CERTIFIED (71/71, terminal gate
# passing, true ±4.5 input set) in globally normalized coordinates, backward-only,
# VOLUME-TUNED. Funnel-size ladder, measured (V = π·√det Q, physical units;
# PR benchmark = Florentin's best sweep config, `trajectory_std` scaling):
#
#   config                          V0 (entry)   Vmin       Vmed      coverage
#   PR best (:logdet+:max_volume)   0.001605     0.000213   0.000767  —
#   ours :maximin,:first_consistent 0.000209     0.000209   0.000267  4021
#   ours :logdet+:max_volume        0.2124       0.01254    0.03673   5.48   ← this file
#
# Normalized coordinates + the volume objective + the box line-search give entry
# funnels 132× the PR's best, with no collapse (Vmin healthy) and the initial-set
# coverage margin down from 4021 to 5.5 — gap D nearly closed as a byproduct.
# Earlier ablations: per-step scaling ⟹ tiny funnels; no scaling ⟹ infeasible at
# k≈63; :maximin equalizes the chain (V0 = Vmin) — the safety choice, not the
# size choice. Remaining: coverage margin 5.5 > 1 (entry semi-axes ~2.3× short of
# the full ±10°×±0.5 initial set).

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
include(joinpath(PROBLEMS, "Pendulum", "simple_pendulum.jl"))

params = SimplePendulum.Params(; l = 1.0, g = 9.81)
problem = SimplePendulum.optimal_control_problem(;
    params = params,
    objective = "reachability_up_medium_power_no_obstacle",
)
Δt = 0.1

periodic_dims = SVector(1)
periods = SVector(2π)
wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = SVector(-π))

# ------------------------------------------------------------
# 1) Abstraction seed
# ------------------------------------------------------------

println("— seed (uniform grid abstraction) —")
seed_time = @elapsed begin
    opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("h"), SVector(3.0 * π / 180.0, 0.1))
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(0.3)),
    )
    MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), Δt)
    MOI.set(opt, MOI.RawOptimizerAttribute("approx_mode"), AB.UniformGridAbstraction.GROWTH)
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("jacobian_bound"),
        SimplePendulum.jacobian_bound(params),
    )
    MOI.set(opt, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
    MOI.set(opt, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
    MOI.set(opt, MOI.RawOptimizerAttribute("periodic_periods"), periods)
    MOI.set(opt, MOI.RawOptimizerAttribute("periodic_start"), SVector(-π))
    MOI.set(opt, MOI.RawOptimizerAttribute("early_stop"), true)
    MOI.set(opt, MOI.Silent(), true)

    seed_gen = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
        opt;
        initial_state = SVector(0.0, 0.0),
        concrete = false,
        nstep = 100,
    )
    AB.set_problem!(seed_gen, problem)
    AB.generate!(seed_gen)
    global seed_traj = AB.get_trajectory(seed_gen)
end
println(
    "  $(round(seed_time; digits = 1)) s, seed states: ",
    seed_traj === nothing ? "NONE" : length(ST.states(seed_traj)),
)

# ------------------------------------------------------------
# 2) MPPI refinement (modern stack)
# ------------------------------------------------------------

println("— MPPI —")
base = PR.discretize_problem(problem, Δt; num_substeps = 5)

# The target box [π−15°, π+15°] straddles the ±π seam: wrapped states near
# upright land on the −π side and would miss its upper half. Split the target
# into its in-period copies so membership (success, truncation, cost) works.
T_split = UT.set_in_period(problem.target_set, periodic_dims, periods, SVector(-π))
discrete_problem = PR.OptimalControlProblem(
    base.system,
    base.initial_set,
    T_split,
    base.state_cost,
    base.transition_cost,
    base.time,
    base.safe_set,
)

u_max = 4.5
# Input reserve is the main lever on funnel size: the feedback image κ(E_k) must
# fit in ±u_max, so whatever the plan does not use, the certificate can.
u_plan = 0.6 * u_max

cost = AB.CompositeCost(
    AB.ReachObjectiveCost(T_split; wrap = wrap),
    AB.InputEffortCost(0.001),
    AB.InputSmoothnessCost(; w_du = 0.05, w_ddu = 0.01),
    AB.DomainPenaltyCost(problem.system.X; wrap = wrap),
)

# CEM update (campaign C1's winner for directed one-shot tasks): elite refit
# handles the multimodal energy-pumping landscape better than softmin averaging.
mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = seed_traj,
    nstep = 90,
    nsamples = 1000,
    niter = 40,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.8)),
    project_input = u -> SVector(clamp(u[1], -u_plan, u_plan)),
    cost = cost,
    wrap_state = (problem, x) -> wrap(x),
    update_rule = :cem,
    elite_frac = 0.05,
    antithetic = true,
)
# x-frame lift: periodic unwrap, shifted so the endpoint lands in the target's
# θ-range.
lift = function (traj)
    lifted = ST.unwrap_trajectory(traj, (1,), (2π,))
    θN = collect(ST.states(lifted))[end][1]
    shift = 2π * round((θN - π) / (2π))
    shift == 0.0 && return lifted
    xs = [SVector(x[1] - shift, x[2]) for x in ST.states(lifted)]
    return ST.Trajectory(xs; inputs = collect(ST.inputs(lifted)))
end

# ------------------------------------------------------------
# 4) Symbolic provider on the true input set
# ------------------------------------------------------------

Symbolics.@variables θ ω τ w1 w2 T
xsym = [θ, ω]
usym = [τ]
wsym = [w1, w2]
f_cont(xloc, uloc) = collect(SimplePendulum.dynamic(params)(xloc, uloc))
f_disc = ST.runge_kutta4(f_cont, xsym, usym, T, 1)
fsymbolic = Symbolics.substitute([f_disc[1] + w1, f_disc[2] + w2], Dict(T => Δt))
Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))
provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    xsym,
    usym,
    wsym,
    [0.0, 0.0],
    ST.format_input_set(problem.system.U),   # the PLANT's input set
    ST.format_noise_set(Wset),
)

sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

# ------------------------------------------------------------
# 5) Globally normalized coordinates (plan.md §4.3): certify in z = x ./ t.
# The scaled dynamics f_z(z,u) = f_x(t.*z, u) ./ t are built symbolically once,
# so the Hessian bounds are EXACT in the working frame — the conditioning benefit
# of the per-step state_scaling without its ~1/σ_min(T)² remainder tax (measured:
# scaling off ⇒ infeasible at k≈63; per-step scaling on ⇒ tiny funnels).
# ------------------------------------------------------------

t = [0.85 * 15.0 * π / 180.0, 0.25]
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{2}(LazySets.low(H) ./ t),
    high = SVector{2}(LazySets.high(H) ./ t),
)
zproblem = PR.OptimalControlProblem(
    base.system,
    zbox(problem.initial_set),
    zbox(problem.target_set),
    nothing,
    nothing,
    base.time,
    nothing,
)

ztraj(traj) = ST.Trajectory(
    [SVector{2}(collect(x) ./ t) for x in ST.states(traj)];
    inputs = collect(ST.inputs(traj)),
)
prepare = traj -> ztraj(lift(traj))

# The unwrapped cover has no θ seam and the objective has no obstacle; θ spans the
# full period, so the domain gate would only re-check |ω| ≤ 7 while flagging the
# lifted θ — disable it here and rely on the generous ω margin. Box radii are in
# z-units (x-units ./ t).
# The :max_volume line search (Florentin's contribution) re-solves each step over
# several box scales and keeps the biggest certified ellipsoid — the volume
# analogue of :first_consistent, ~4× the SDPs per step.
adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    true,
    [0.05, 0.10] ./ t,
    [0.005, 0.005] ./ t,
    [2.5, 3.5] ./ t,
    [0.25],
    [0.01],
    [4.5],
    1.5,
    1.05,
    30,
    1e-8,
    false,
    [0.75, 1.0, 1.5, 2.0],
    :max_volume,
    false,
)
back_opts = EB.ChainOptions(;
    maxδx = 12.0,
    maxδu = 3.0,
    λ = 0.001,
    terminal_shape = nothing,            # default: inscribed ellipsoid of the target
    terminal_shrink = 0.95,
    state_scaling = nothing,             # exact: the dynamics are already normalized
    linearization_δx = [0.05, 0.10] ./ t,
    linearization_δu = [1.0],
    adaptive_boxes = adaptive_opts,
    # :logdet maximizes true volume — the robustness objective; the collapse risk
    # it carries (vs :maximin) is monitored by the funnel-area stats below.
    objective = :logdet,
    check_state_domain = false,
)

println("— generate ⇄ certify loop (retry ladder, up to 5 rounds) —")
bw = EB.BackwardCertifier(zprovider, sdp, back_opts)
driver = AB.TrajectoryCertificationOptimizer.Optimizer(mppi, bw)
MOI.set(driver, MOI.RawOptimizerAttribute("concrete_problem"), discrete_problem)
MOI.set(driver, MOI.RawOptimizerAttribute("certifier_problem"), zproblem)
MOI.set(driver, MOI.RawOptimizerAttribute("max_rounds"), 5)
MOI.set(driver, MOI.RawOptimizerAttribute("prepare_trajectory"), prepare)
loop_time = @elapsed MOI.optimize!(driver)

rounds = MOI.get(driver, MOI.RawOptimizerAttribute("rounds"))
loop_success = MOI.get(driver, MOI.RawOptimizerAttribute("success"))
traj = MOI.get(driver, MOI.RawOptimizerAttribute("trajectory"))
bres = EB.get_result(bw)
bres === nothing &&
    error("the generator failed in every round — no certification was attempted")
d = AB.MPPITrajectoryGenerator.get_diagnostics(mppi)
println("  $(round(loop_time; digits = 1)) s, rounds = $rounds, certified = $loop_success")
println(
    "  last round: mppi steps = $(length(ST.inputs(traj))), ",
    "backward failed_k = $(bres.failed_k), terminal ⊆ target: ",
    "$(bres.terminal_contained), initial coverage: $(bres.initial_coverage)",
)
if !bres.success && !isempty(bres.steps)
    last_step = bres.steps[1]
    println(
        "  first failing record: k = $(last_step.k), ",
        get(
            last_step.summary,
            :gate_reason,
            get(last_step.summary, :adaptive_box_status, "?"),
        ),
    )
end
lifted_x = lift(traj)

# Backward-only: the certificate needs no forward chain — the funnel controller
# below is the complete product. (`κ_z(z) = K_z·z + b` maps back to the physical
# frame as `κ_x(x) = K_z·D⁻¹·x + b`.)
# Funnel sizes in the physical frame, comparable to the PR's sweep report
# (V = π·√det(Q), the ellipse area; his best: V0 = 1.6e-3, Vmin = 2.1e-4).
if bres.success
    areas = [
        π * sqrt(LA.det(D * Matrix(LazySets.shape_matrix(E)) * D)) for
        E in bres.lmi_data.ellipsoids
    ]
    sorted_areas = sort(areas)
    println(
        "  funnel areas: V0 = $(round(areas[1]; sigdigits = 4)), ",
        "Vmin = $(round(sorted_areas[1]; sigdigits = 4)), ",
        "Vmed = $(round(sorted_areas[div(end, 2)]; sigdigits = 4)), ",
        "Vmax = $(round(sorted_areas[end]; sigdigits = 4))",
    )
end

if loop_success
    zctrl = AB.get_controller(bw)
    ctrl = ST.FunnelController(
        [MS.AffineMap(Matrix(κ.A) / Matrix(D), collect(κ.c)) for κ in zctrl.kappas],
        [
            LazySets.Ellipsoid(
                t .* collect(LazySets.center(E)),
                Matrix(D * Matrix(LazySets.shape_matrix(E)) * D),
            ) for E in zctrl.ellipsoids
        ],
    )
    println(
        "— certified controller (physical frame): FunnelController with ",
        "$(length(ctrl.kappas)) steps —",
    )
end

# ------------------------------------------------------------
# 6) Plot: trajectory + certified sets in the lifted (θ, ω) plane
# ------------------------------------------------------------

fig = plot(;
    xlabel = "θ  [rad]",
    ylabel = "ω  [rad/s]",
    title = "Certified pendulum swing-up",
    legend = :topleft,
    size = (900, 600),
)

plot!(fig, problem.initial_set; color = :gray, alpha = 0.5, label = "initial set")
plot!(fig, problem.target_set; color = :green, alpha = 0.35, label = "target set")

# De-normalize the certified funnel back to the physical frame for the plot.
denorm(E) = LazySets.Ellipsoid(
    t .* collect(LazySets.center(E)),
    Matrix(D * Matrix(LazySets.shape_matrix(E)) * D),
)

funnel = [denorm(E) for E in bres.lmi_data.ellipsoids]
for (i, E) in enumerate(funnel)
    plot!(
        fig,
        E;
        color = :steelblue,
        alpha = 0.25,
        linewidth = 1.2,
        linecolor = :steelblue,
        label = i == 1 ? "backward funnel" : "",
    )
end

xs_plot = collect(ST.states(lifted_x))
plot!(
    fig,
    [x[1] for x in xs_plot],
    [x[2] for x in xs_plot];
    color = :black,
    linewidth = 2,
    marker = :circle,
    markersize = 2,
    label = "nominal trajectory",
)

# Zoom panel: the funnel up close on the mid-swing segment, where the ellipsoids
# are actually visible at scale.
mid = xs_plot[div(length(xs_plot), 2)]
figz = plot(;
    xlabel = "θ  [rad]",
    ylabel = "ω  [rad/s]",
    title = "zoom: mid-swing ellipsoids",
    legend = false,
    xlims = (mid[1] - 1.0, mid[1] + 1.0),
    ylims = (mid[2] - 1.6, mid[2] + 1.6),
)
for E in funnel
    plot!(figz, E; color = :steelblue, alpha = 0.3, linewidth = 1.5, linecolor = :navy)
end
plot!(
    figz,
    [x[1] for x in xs_plot],
    [x[2] for x in xs_plot];
    color = :black,
    linewidth = 2,
    marker = :circle,
    markersize = 3,
)

final = plot(fig, figz; layout = (1, 2), size = (1500, 620))
plot_path = joinpath(@__DIR__, "demo_pendulum.png")
savefig(final, plot_path)
println("— plot saved: $plot_path —")
