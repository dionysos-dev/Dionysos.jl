# Certified parallel parking (créneau) — the master-thesis benchmark from
# PR #569, on the refactored pipeline.
#
# 4-D tractor-trailer with a realistic rig (L1 = 2.2, L2 = 2.8, Lc = 0.9) on a
# 21 × 5.6 m street: two parked cars along the curb leave a 9 m slot. The rig
# starts in the traffic lane AHEAD of the slot heading along the street and
# reverses trailer-first through the classic S-curve back into the slot,
# ending parallel-parked (θ ≈ 0). Unlike the reverse-bay demo, θ never
# approaches ±π — no periodic cover, no unwrap.
#
#     julia --project=test research/TrajectoryCertificationOptimizer/demo_vehicle_creneau.jl
#
# Pipeline: scripted S-curve seed scan (straight reverse → arc → counter-arc →
# feedback settle, domain-clean scripts only) → CEM-MPPI → backward :maximin
# certification in normalized coordinates with the state-domain gate ACTIVE
# against the parked cars → capture-phase fallback → closed-loop
# falsification.
#
# Status (2026-08-14, Clarabel, fixed rng): FULL-CHAIN certificate, round 1 —
# the entire 70-step maneuver certifies end-to-end (failed_k = nothing,
# terminal ⊆ target, wall-aware domain gates green), no capture fallback
# needed; endpoint (11.30, 1.05, 0.00, 0.00) = the slot center exactly
# (box-centered terminal engaged); closed-loop 50/50. The funnels are thin
# tubes (V_entry ≈ 2e-11): reversing an articulated rig is open-loop unstable,
# so certified tubes pay the same authority price as demo_vehicle — and
# :logdet pancake-collapses in 4-D there (measured), which is why this demo
# ships :maximin.
#
# Body-clearance lessons (each measured via the rig-body sweep probe):
# - The POINT domain hides the body: the trailer ends ≈ 3.9 m behind the
#   state, so the rear car carries a TRAILER-SHADOW extension (east face
#   x = 10.5) and the target's faces are set by the parked BODY (left 10.6 for
#   the trailer, right 12.0 for the 2.75 m tractor nose).
# - The PR's 8 m slot leaves 1.3 m of spare for the 6.7 m rig — parkable only
#   by multi-shunt maneuvers; the benchmark ships a 9 m slot (2.3 m spare),
#   the tightest a one-shot S enters cleanly.
# - The terminal pull is a SEARCH metric: with the true ±5.5° angular radii
#   the seed's residual θ outprices a 1.5 m altitude error and CEM settles
#   into a hover ABOVE the slot; ±15° radii keep the basin pointing down.

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
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
include(joinpath(PROBLEMS, "ArticulatedVehicle", "articulated_vehicle.jl"))
const AV = ArticulatedVehicle

params = AV.parallel_parking_params()
# The benchmark's collision model is the tractor rear axle POINT; the drawn
# bodies are wider AND the trailer ends Lc + L2 (+ wheels) ≈ 3.9 m BEHIND the
# point. Two corrections make the certified point trajectory keep the drawn
# rig clear of the TRUE cars (plots show the true cars):
# - a lateral body-clearance margin (trailer half-width 0.44 m + slack), and
# - the rear car's TRAILER SHADOW: in the low-y slot band a tractor west of
#   x = 10.5 parks the trailer INTO the rear car (measured: the un-shadowed
#   optimum dives to x ≈ 10.3 mid-settle, trailer kissing the car), so the
#   rear obstacle extends east to 10.5 for the domain.
obstacle_margin = 0.5
domain_obstacles = [
    # Rear car + TRAILER shadow: the trailer ends Lc + L2 (+ wheels) ≈ 3.9 m
    # BEHIND the tractor point, so in the low-y slot band the point must stay
    # east of x = 10.5 for the drawn trailer to clear the TRUE rear car
    # (measured: without the shadow the optimum dives to x ≈ 10.3, trailer
    # kissing the car). The east face is a LONGITUDINAL bound on the point and
    # takes no lateral margin. The tractor NOSE (2.75 m ahead) needs no shadow
    # in the 9 m slot: the target's right face (12.0) already keeps it ≥ 0.55 m
    # clear of the front car.
    LazySets.Hyperrectangle(;
        low = SVector(-obstacle_margin, -obstacle_margin),
        high = SVector(10.5, 2.1 + obstacle_margin),
    ),
    AV._inflate(
        LazySets.Hyperrectangle(; low = SVector(15.3, 0.0), high = SVector(21.0, 2.1)),
        obstacle_margin,
    ),
]
problem = AV.parallel_parking_problem(; params = params, obstacles2d = domain_obstacles)

Δt = 0.2
base = PR.discretize_problem(problem, Δt; num_substeps = 4)
f = MS.mapping(base.system)

# ------------------------------------------------------------
# 1) Seed: S-curve reverse script scan. Phases: straight reverse (n0), arc
# (δ1, n1 — nose right, rear swings toward the curb), counter-arc (−δ2, n2),
# then a feedback settle (n3) with the trailer-backing law δ = −2ϕ − θ that
# straightens the rig inside the slot.
# ------------------------------------------------------------

x0 = SVector(19.3, 4.0, 0.0, 0.0)
x_goal = SVector(11.3, 1.05, 0.0, 0.0)         # center of the slot target

function rollout_script(n0, δ1, n1, δ2, n2, n3, n4)
    xs = [x0]
    us = SVector{2, Float64}[]
    step!(u) = (push!(us, u); push!(xs, f(xs[end], u)))
    for _ in 1:n0
        step!(SVector(-1.0, 0.0))
    end
    for _ in 1:n1
        step!(SVector(-1.0, δ1))
    end
    for _ in 1:n2
        step!(SVector(-1.0, -δ2))
    end
    for _ in 1:n3
        x = xs[end]
        step!(SVector(-1.0, clamp(-2.0 * x[4] - x[3], -0.85, 0.85)))
    end
    # Forward correction — the human move: back up to the shadow limit, then
    # pull forward to center in the slot (forward driving is the ϕ-stable
    # direction, so a simple heading feedback settles the rig).
    for _ in 1:n4
        x = xs[end]
        step!(SVector(1.0, clamp(-1.5 * x[3], -0.5, 0.5)))
    end
    return xs, us
end

score(x) = (x[1] - x_goal[1])^2 + (x[2] - x_goal[2])^2 + 4 * x[3]^2 + 4 * x[4]^2

# Asymmetric S (δ1 ≠ δ2): under the trailer-shadow domain the settle cannot
# dive west to fix the alignment, so the arcs themselves must place y and θ.
best = (Inf, 0, 0.4, 10, 0.4, 10, 8, 0)
for n0 in 0:3:12,
    δ1 in 0.3:0.1:0.8,
    n1 in 6:2:20,
    δ2 in 0.3:0.1:0.8,
    n2 in 6:2:20,
    n3 in 0:6:24,
    n4 in 0:6:18

    xs, _ = rollout_script(n0, δ1, n1, δ2, n2, n3, n4)
    all(x -> collect(x) ∈ problem.system.X, xs) || continue
    s = score(xs[end])
    s < best[1] && (global best = (s, n0, δ1, n1, δ2, n2, n3, n4))
end
isfinite(best[1]) || error("seed scan: no domain-clean script found")
_, n0, δ1, n1, δ2, n2, n3, n4 = best
println(
    "— seed scan: n0 = $n0, arc δ = $δ1 × $n1, counter-arc × $n2, settle × $n3, " *
    "forward × $n4, endpoint score = $(round(best[1]; digits = 3)) —",
)

nstep = n0 + n1 + n2 + n3 + n4 + 12            # slack for CEM to reshape
seed_traj = begin
    xs, us = rollout_script(n0, δ1, n1, δ2, n2, n3, n4)
    while length(us) < nstep                   # v = 0 freezes the rig in place
        push!(us, SVector(0.0, 0.0))
        push!(xs, f(xs[end], us[end]))
    end
    ST.Trajectory(xs; inputs = us)
end

# ------------------------------------------------------------
# 2) CEM-MPPI
# ------------------------------------------------------------

println("— MPPI —")
u_plan = SVector(2.0, 0.6)                     # reserve vs the plant's [±5, ±0.96]

cost = AB.CompositeCost(
    AB.ReachObjectiveCost(problem.target_set),
    # The pull is a SEARCH metric, not the target box: with tight angular radii
    # the seed's residual θ error outprices a 1.5 m altitude error and CEM
    # abandons the slot for a hover above it (measured). Wide angles keep the
    # basin pointing DOWN; the reach cost still enforces the true box.
    AB.TerminalPullCost(
        collect(x_goal),
        [0.63, 0.5, deg2rad(15.0), deg2rad(15.0)];
        w = 1e4,
    ),
    AB.InputEffortCost(0.01),
    AB.InputSmoothnessCost(; w_du = 0.5, w_ddu = 0.1),
    # Dominates the terminal pull so domain adherence is never traded away —
    # certification dies wherever the nominal leaves the domain.
    AB.DomainPenaltyCost(problem.system.X; w = 1e5),
)

mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = seed_traj,
    nstep = nstep,
    nsamples = 800,
    niter = 180,
    anneal = 0.99,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.12, 0.08)),
    project_input = u ->
        SVector(clamp(u[1], -u_plan[1], u_plan[1]), clamp(u[2], -u_plan[2], u_plan[2])),
    cost = cost,
    update_rule = :cem,
    elite_frac = 0.05,
    antithetic = true,
    stop_on_success = false,
)

# ------------------------------------------------------------
# 3) Normalized symbolic provider
# ------------------------------------------------------------

Symbolics.@variables xv[1:4] uv[1:2] wv[1:4] T
xsym = collect(xv)
usym = collect(uv)
wsym = collect(wv)
f_cont(xloc, uloc) = collect(AV.dynamic(params)(xloc, uloc))
f_disc = ST.runge_kutta4(f_cont, xsym, usym, T, 1)
fsymbolic = Symbolics.substitute([f_disc[i] + wsym[i] for i in 1:4], Dict(T => Δt))
Wset = LazySets.Hyperrectangle(;
    low = zero(SVector{4, Float64}),
    high = zero(SVector{4, Float64}),
)
provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    xsym,
    usym,
    wsym,
    zeros(4),
    ST.format_input_set(problem.system.U),     # the plant's true input set
    ST.format_noise_set(Wset),
)

t = [3.0, 1.5, 0.4, 0.4]                       # characteristic scales (m, m, rad, rad)
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{4}(LazySets.low(H) ./ t),
    high = SVector{4}(LazySets.high(H) ./ t),
)

street_box = LazySets.Hyperrectangle(;
    low = SVector(0.0, 0.0, -pi, -deg2rad(65.0)),
    high = SVector(21.0, 5.6, pi, deg2rad(65.0)),
)
zX = UT.set_minus(
    zbox(street_box),
    UT.set_union([
        zbox(AV.extrude_xy_obstacle_to_4d(ob, street_box)) for ob in domain_obstacles
    ]),
)
zsys = MS.ConstrainedBlackBoxControlContinuousSystem(
    AV.dynamic(params),
    4,
    2,
    zX,
    problem.system.U,
)
zproblem = PR.OptimalControlProblem(
    zsys,
    zbox(problem.initial_set),
    zbox(problem.target_set),
    nothing,
    nothing,
    base.time,
    nothing,
)

prepare = traj -> ST.normalize_trajectory(traj, t)

# ------------------------------------------------------------
# 4) Backward certification, domain gate ON (parked cars included)
# ------------------------------------------------------------

sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

adaptive_opts = EB.AdaptiveLinearizationBoxOptions(;
    ΔX_initial = [0.05, 0.05, 0.05, 0.05] ./ t,
    ΔX_min = [0.005, 0.005, 0.005, 0.005] ./ t,
    ΔX_max = [3.0, 3.0, 1.0, 1.0] ./ t,
    ΔU_initial = [0.3, 0.1],
    ΔU_min = [0.03, 0.01],
    ΔU_max = [2.0, 0.6],
)
back_opts = EB.ChainOptions(;
    maxδx = 8.0,
    maxδu = 3.0,
    λ = 0.001,
    terminal_shape = nothing,
    terminal_shrink = 0.9,
    linearization_δx = [0.05, 0.05, 0.05, 0.05] ./ t,
    linearization_δu = [0.3, 0.1],
    adaptive_boxes = adaptive_opts,
    objective = :maximin,
    remainder_model = :vertices,               # n = 4: 16 blocks/step is affordable
    domain_cap = true,
    check_state_domain = true,
)

println("— generate ⇄ certify loop (retry ladder, up to 3 rounds) —")
bw = EB.BackwardCertifier(zprovider, sdp, back_opts)
driver = AB.TrajectoryCertificationOptimizer.Optimizer(mppi, bw)
MOI.set(driver, MOI.RawOptimizerAttribute("concrete_problem"), base)
MOI.set(driver, MOI.RawOptimizerAttribute("certifier_problem"), zproblem)
MOI.set(driver, MOI.RawOptimizerAttribute("max_rounds"), 3)
MOI.set(driver, MOI.RawOptimizerAttribute("prepare_trajectory"), prepare)
loop_time = @elapsed MOI.optimize!(driver)

rounds = MOI.get(driver, MOI.RawOptimizerAttribute("rounds"))
loop_success = MOI.get(driver, MOI.RawOptimizerAttribute("success"))
traj = MOI.get(driver, MOI.RawOptimizerAttribute("trajectory"))
bres = EB.get_result(bw)
bres === nothing &&
    error("the generator failed in every round — no certification was attempted")
println("  $(round(loop_time; digits = 1)) s, rounds = $rounds, certified = $loop_success")
println(
    "  last round: mppi steps = $(length(ST.inputs(traj))), ",
    "backward failed_k = $(bres.failed_k), terminal ⊆ target: ",
    "$(bres.terminal_contained), initial coverage: $(bres.initial_coverage), ",
    "domain gate ran: $(bres.state_domain_checked)",
)
if !bres.success && !isempty(bres.steps)
    first_bad = bres.steps[1]
    println(
        "  first failing record: k = $(first_bad.k), ",
        get(
            first_bad.summary,
            :gate_reason,
            get(first_bad.summary, :adaptive_box_status, "?"),
        ),
    )
end

# Nominal-in-domain check: wherever the nominal leaves the domain, the per-step
# domain cap silently disengages and the gate rejects a posteriori.
let xs = collect(ST.states(traj)), out = findall(x -> collect(x) ∉ problem.system.X, xs)
    if !isempty(out)
        println("  WARNING: nominal outside the domain at state(s) $out:")
        for i in out
            println("    x[$i] = $(round.(collect(xs[i]); digits = 3))")
        end
    end
end
let xe = collect(ST.states(traj))[end]
    r_box = 0.9 .* [0.7, 0.55, deg2rad(6.0), deg2rad(3.0)]
    margin = sum(((collect(xe) .- collect(x_goal)) ./ r_box) .^ 2)
    println(
        "  endpoint = $(round.(collect(xe); digits = 3)); box-centered margin = ",
        "$(round(margin; digits = 3)) (engages iff ≤ 0.64)",
    )
end

# ------------------------------------------------------------
# 4b) Capture-phase fallback: re-certify the certified suffix standalone.
# ------------------------------------------------------------

if !loop_success && bres.failed_k !== nothing && !isempty(bres.lmi_data.ellipsoids)
    ksuf = bres.failed_k + 1
    xs_t = collect(ST.states(traj))
    us_t = collect(ST.inputs(traj))
    suffix = ST.Trajectory(xs_t[ksuf:end]; inputs = us_t[ksuf:end])
    println(
        "— capture-phase certificate: re-certifying the suffix (states $ksuf..$(length(xs_t))) —",
    )
    AB.set_trajectory!(bw, ST.normalize_trajectory(suffix, t))
    suffix_time = @elapsed AB.certify!(bw)
    bres = EB.get_result(bw)
    println(
        "  $(round(suffix_time; digits = 1)) s, success = $(bres.success), ",
        "failed_k = $(bres.failed_k), terminal ⊆ target: $(bres.terminal_contained)",
    )
    global loop_success = bres.success
end

ells = bres.lmi_data.ellipsoids
if !isempty(ells)
    vols = [π^2 / 2 * sqrt(LA.det(D * Matrix(LazySets.shape_matrix(E)) * D)) for E in ells]
    sv = sort(vols)
    println(
        "  funnel volumes (certified segment): V_entry = $(round(vols[1]; sigdigits = 4)), ",
        "Vmin = $(round(sv[1]; sigdigits = 4)), Vmed = $(round(sv[div(end, 2)]; sigdigits = 4)), ",
        "Vmax = $(round(sv[end]; sigdigits = 4))",
    )
end

loop_success && println(
    "— certified controller: FunnelController with $(length(AB.get_controller(bw).kappas)) steps —",
)

# ------------------------------------------------------------
# 5) Closed-loop falsification + inflation stress (α = 1 is the certificate;
# α > 1 measures the empirical success basin beyond it, failures by mode).
# ------------------------------------------------------------

if !isempty(ells) && length(ells) > 1
    i_val = min(something(findfirst(v -> v >= 1e-6, vols), 1), length(ells) - 1)
    f_z(z, u) = SVector{4}(f(SVector{4}(t .* z), u) ./ t)
    rows = EB.inflation_stress(
        f_z,
        bres.lmi_data.kappas[i_val:end],
        bres.lmi_data.ellipsoids[i_val:end],
        zbox(problem.target_set);
        alphas = [1.0, 1.5, 2.0, 3.0],
        n_samples = 100,
        rng = Random.MersenneTwister(7),
        input_set = problem.system.U,
        project_input = u ->
            SVector(clamp(u[1], -5.0, 5.0), clamp(u[2], -0.959931, 0.959931)),
        domain = zX,
    )
    println(
        "  closed-loop + inflation stress from the ",
        "$(round((length(ells) - i_val) * Δt; digits = 2))-s-out region ",
        "(V = $(round(vols[i_val]; sigdigits = 3))):",
    )
    EB.print_inflation_stress(rows)
end

# ------------------------------------------------------------
# 6) Plot: street with the parked cars, trajectory, funnel shadows + slot zoom
# ------------------------------------------------------------

fig = plot(;
    xlabel = "x  [m]",
    ylabel = "y  [m]",
    title = "Certified parallel parking (tractor-trailer)",
    legend = :topleft,
    size = (1000, 480),
    aspect_ratio = :equal,
    xlims = (-0.3, 21.3),
    ylims = (-0.3, 6.0),
)

box2(H) = LazySets.Hyperrectangle(;
    low = SVector(LazySets.low(H, 1), LazySets.low(H, 2)),
    high = SVector(LazySets.high(H, 1), LazySets.high(H, 2)),
)
AV.plot_xy_obstacles!(fig, AV.parallel_parking_obstacles())
AV.plot_xy_obstacles!(
    fig,
    domain_obstacles;
    alpha = 0.0,
    linestyle = :dash,
    linecolor = :gray,
)
plot!(fig, box2(problem.initial_set); color = :gray, alpha = 0.5, label = "initial set")
plot!(fig, box2(problem.target_set); color = :green, alpha = 0.35, label = "target (slot)")

funnel_shadow(E) = begin
    Q = D * Matrix(LazySets.shape_matrix(E)) * D
    c = t .* collect(LazySets.center(E))
    LazySets.Ellipsoid(c[1:2], LA.Symmetric(Q[1:2, 1:2]) |> Matrix)
end
shadows = [funnel_shadow(E) for E in ells]
for (i, E) in enumerate(shadows)
    plot!(
        fig,
        E;
        color = :steelblue,
        alpha = 0.25,
        linewidth = 1.0,
        linecolor = :steelblue,
        label = i == 1 ? "funnel (xy shadow)" : "",
    )
end

xs_plot = collect(ST.states(traj))
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

figz = plot(;
    xlabel = "x  [m]",
    ylabel = "y  [m]",
    title = "zoom: settling into the slot",
    legend = false,
    aspect_ratio = :equal,
    xlims = (8.0, 16.0),
    ylims = (0.0, 3.2),
)
AV.plot_xy_obstacles!(figz, AV.parallel_parking_obstacles())
plot!(figz, box2(problem.target_set); color = :green, alpha = 0.35)
for E in shadows
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

final = plot(fig, figz; layout = (1, 2), size = (1500, 520))
plot_path = joinpath(@__DIR__, "demo_vehicle_creneau.png")
savefig(final, plot_path)
println("— plot saved: $plot_path —")

# ------------------------------------------------------------
# 7) Dashboard animation: rig in the street + xy and input panels
# ------------------------------------------------------------

av_plot = AV.system_plot!(;
    params = params,
    obstacles2d = AV.parallel_parking_obstacles(),
    xlims = (-0.5, 21.5),
    ylims = (-2.0, 7.0),
)
dash_plot! =
    (f, xk, uk) -> begin
        plot!(f, box2(problem.initial_set); color = :gray, alpha = 0.4)
        plot!(f, box2(problem.target_set); color = :green, alpha = 0.4)
        av_plot(f, xk, uk)
    end
dash_path = DI.animate_trajectory_dashboard(
    dash_plot!,
    traj;
    xdims = (1, 2),
    udims = (1, 2),
    Δt = Δt,
    fps = 12,
    filename = joinpath(@__DIR__, "demo_vehicle_creneau_dashboard.gif"),
    title = "Certified parallel parking",
)
println("— dashboard saved: $dash_path —")
