# Certified reverse parking — the "Permis" benchmark from PR #569, on the
# refactored pipeline.
#
# 4-D tractor-trailer (x, y, θ, ϕ), compact rig (L1 = 1, L2 = 1, Lc = 0.5).
# Two walls leave a dead-end bay (corridor y ∈ [4.7, 6.0] for x ∈ [4, 10])
# opening west. The vehicle starts near the origin heading east, must U-turn in
# the open area and BACK the trailer into the bay, ending θ ≈ π (facing back
# out) with the rig straight — reversing trailer-first is open-loop
# jackknife-unstable, which is exactly what the certified funnel feedback
# stabilizes.
#
#     julia --project=test research/TrajectoryCertificationOptimizer/demo_vehicle_parking.jl
#
# Pipeline: scripted two-phase seed scan (forward U-turn arc + feedback-
# stabilized reverse) → CEM-MPPI on the θ-extended cover → periodic unwrap +
# shift → backward :maximin certification in normalized coordinates with the
# state-domain gate ACTIVE against the walls → capture-phase fallback →
# closed-loop falsification.
#
# Status (2026-08-13, Clarabel, fixed rng): generation SOLVED — 50 steps, a
# 184° U-turn (δ = 0.3) then a 24-step stabilized reverse into the bay,
# endpoint (9.11, 5.51, π, 0.0), terminal ⊆ target, nominal fully in-domain.
# The full backward chain certifies 50 → 15 and dies at k = 14 with
# `lmi_infeasible_at_max_box`: funnel volume bleeds ~2×/step backward through
# the reverse (the certified direction fights the jackknife instability), so
# by the benign forward arc the target ellipsoids are needles the remainder
# swamps — the same wall class as demo_vehicle's mid-turn. SHIPPED result:
# the CAPTURE-PHASE CERTIFICATE — the suffix (states 15..51, i.e. the ENTIRE
# jackknife-unstable reverse phase) re-certified standalone as a complete
# chain (36-step FunnelController, terminal + wall-aware domain gates green),
# validated closed-loop 50/50.
#
# Measured traps (each cost a debugging round):
# - Open-loop reverse jackknifes: with δ = 0 at v = −2, ϕ̇ ≈ 2 sin ϕ — a
#   constant-input reverse script ends at |ϕ| ≈ π. Even the SEED needs the
#   backing feedback δ = −2ϕ − (θ − π) (linearly stable: trace −2, det 4).
# - The seed scan must reject domain-dirty scripts: an over-rotated arc (217°)
#   scores a great ENDPOINT — the reverse feedback tracks back — while its
#   U-turn swings 0.7 m west of the x = −1 face. Wherever the nominal leaves
#   the domain, `domain_cap` silently disengages and the funnel dies at the
#   state-domain gate (k = 28, axis-1 support excess 0.165 — found by the
#   gate's per-axis violation report).
# - The θ cover must be ±2.5π, not ±2π: a right-turning solution (θ → −π)
#   shifted +2π onto the target's +π branch parks its early states exactly on
#   a ±2π cover face.

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

params = AV.parking_params()
problem = AV.parking_problem(; params = params)
# θ-extended cover: MPPI rollouts and the certified chain live on the unwrapped
# θ axis, so every domain check — the CEM penalty and the certifier's
# state-domain gate — runs against the cover with the walls extruded over the
# full θ range. ±2.5π, NOT ±2π: a right-turning solution (θ → −π) gets shifted
# by +2π onto the target's +π branch, which parks its early states (θ ≈ 0 with
# wobble) exactly on a ±2π face — there the nominal leaves the box, the domain
# cap silently disengages, and the funnel dies at the gate (measured).
const θ_cov = 2.5pi
cover = AV.parking_problem(; params = params, theta_lims = (-θ_cov, θ_cov))

Δt = 0.2
base = PR.discretize_problem(problem, Δt; num_substeps = 4)
base_cover = PR.discretize_problem(cover, Δt; num_substeps = 4)
f = MS.mapping(base.system)

# Both unwrapped branches of the target (θ ≈ π and θ ≈ -π): generation success,
# truncation, and the reach cost must accept either turn direction.
T_plus = problem.target_set
T_minus = LazySets.Hyperrectangle(;
    low = SVector(9.0, 5.0, -pi - deg2rad(5.0), -deg2rad(5.0)),
    high = SVector(10.0, 6.0, -pi + deg2rad(5.0), deg2rad(5.0)),
)
T_union = UT.set_union([T_plus, T_minus])
gen_problem = PR.OptimalControlProblem(
    base_cover.system,
    base.initial_set,
    T_union,
    base.state_cost,
    base.transition_cost,
    base.time,
    nothing,
)

# ------------------------------------------------------------
# 1) Seed: two-phase script scan (forward U-turn arc, then straight reverse)
# ------------------------------------------------------------

x0 = SVector(0.0, 0.0, 0.0, 0.0)
x_goal = SVector(9.5, 5.5, pi, 0.0)

# Phase 2 reverses under the classic trailer-backing feedback
# δ = −k_ϕ·ϕ − k_θ·(θ−π): open-loop reverse jackknifes (ϕ̇ ≈ 2 sin ϕ at v = −2,
# unstable at 0 — measured: a constant-δ script ends at |ϕ| ≈ π), while this
# law is linearly stable (trace −2, det 4). The recorded inputs still make a
# plain open-loop seed for CEM.
function rollout_script(δ_arc, n1, n2)
    xs = [x0]
    us = SVector{2, Float64}[]
    for _ in 1:n1
        push!(us, SVector(2.0, δ_arc))
        push!(xs, f(xs[end], us[end]))
    end
    for _ in 1:n2
        x = xs[end]
        δ = clamp(-2.0 * x[4] - 1.0 * (x[3] - pi), -0.45, 0.45)
        push!(us, SVector(-2.0, δ))
        push!(xs, f(xs[end], us[end]))
    end
    return xs, us
end

score(x) = begin
    dθ = mod(x[3] - x_goal[3] + pi, 2pi) - pi
    (x[1] - x_goal[1])^2 + (x[2] - x_goal[2])^2 + 4 * dθ^2 + 4 * x[4]^2
end

# Domain adherence is part of the scan: an over-rotated arc (e.g. δ = 0.35 for
# 26 steps = 217°) scores a great ENDPOINT — the reverse feedback tracks back —
# while its U-turn swings 0.7 m west of the x = −1 face; certification then
# dies where the nominal leaves the domain (measured). Scripts are scored on
# the cover (raw θ is unwrapped).
best = (Inf, 0.35, 22, 24)
for δ_arc in 0.25:0.05:0.5, n1 in 16:2:36, n2 in 16:2:36
    xs, _ = rollout_script(δ_arc, n1, n2)
    all(x -> collect(x) ∈ cover.system.X, xs) || continue
    s = score(xs[end])
    s < best[1] && (global best = (s, δ_arc, n1, n2))
end
isfinite(best[1]) || error("seed scan: no domain-clean script found")
_, δ_arc, n1, n2 = best
println(
    "— seed scan: δ_arc = $δ_arc, n1 = $n1 (U-turn), n2 = $n2 (reverse), " *
    "endpoint score = $(round(best[1]; digits = 3)) —",
)

nstep = n1 + n2 + 14                           # slack for CEM to reshape both phases
seed_traj = begin
    xs, us = rollout_script(δ_arc, n1, n2)
    while length(us) < nstep                   # v = 0 freezes the rig in place
        push!(us, SVector(0.0, 0.0))
        push!(xs, f(xs[end], us[end]))
    end
    ST.Trajectory(xs; inputs = us)
end

# ------------------------------------------------------------
# 2) CEM-MPPI on the cover
# ------------------------------------------------------------

println("— MPPI —")
u_plan = SVector(2.5, 0.45)                    # reserve vs the plant's [±5, ±π/4]

cost = AB.CompositeCost(
    AB.ReachObjectiveCost(T_union),
    # Periodic pull to the bay pose: distance-normalized by the target radii,
    # θ measured on the circle so either turn direction is rewarded.
    AB.TerminalPullCost(
        collect(x_goal),
        [0.45, 0.45, deg2rad(4.5), deg2rad(4.5)];
        w = 3e4,
        periods = [nothing, nothing, 2pi, nothing],
    ),
    AB.InputEffortCost(0.01),
    AB.InputSmoothnessCost(; w_du = 0.5, w_ddu = 0.1),
    # Dominates the terminal pull so domain adherence is never traded away —
    # certification dies wherever the nominal leaves the domain (see header).
    AB.DomainPenaltyCost(cover.system.X; w = 1e5),
)

mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = seed_traj,
    nstep = nstep,
    nsamples = 800,
    niter = 120,
    anneal = 0.99,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.15, 0.06)),
    project_input = u ->
        SVector(clamp(u[1], -u_plan[1], u_plan[1]), clamp(u[2], -u_plan[2], u_plan[2])),
    cost = cost,
    update_rule = :cem,
    elite_frac = 0.05,
    antithetic = true,
    stop_on_success = false,
)

# Unwrapped lift, shifted so the endpoint lands on the +π branch of the target.
lift = function (traj)
    lifted = ST.unwrap_trajectory(traj, (3,), (2pi,))
    θN = collect(ST.states(lifted))[end][3]
    shift = 2pi * round((θN - pi) / (2pi))
    shift == 0.0 && return lifted
    xs = [SVector(x[1], x[2], x[3] - shift, x[4]) for x in ST.states(lifted)]
    return ST.Trajectory(xs; inputs = collect(ST.inputs(lifted)))
end

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

t = [2.0, 2.0, 0.35, 0.35]                     # characteristic scales (m, m, rad, rad)
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{4}(LazySets.low(H) ./ t),
    high = SVector{4}(LazySets.high(H) ./ t),
)

# z-frame cover domain WITH the walls: scale the cover box and each extruded
# obstacle, subtract — the state-domain gate then proves every funnel ellipsoid
# disjoint from the walls in the frame it was certified in.
cover_box = LazySets.Hyperrectangle(;
    low = SVector(-1.0, -1.0, -θ_cov, -deg2rad(50.0)),
    high = SVector(10.0, 9.0, θ_cov, deg2rad(50.0)),
)
zX = UT.set_minus(
    zbox(cover_box),
    UT.set_union([
        zbox(AV.extrude_xy_obstacle_to_4d(ob, cover_box)) for ob in AV.parking_obstacles()
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
    zbox(T_plus),
    nothing,
    nothing,
    base.time,
    nothing,
)

prepare = traj -> ST.normalize_trajectory(lift(traj), t)

# ------------------------------------------------------------
# 4) Backward certification, domain gate ON (walls included)
# ------------------------------------------------------------

sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

adaptive_opts = EB.AdaptiveLinearizationBoxOptions(;
    ΔX_initial = [0.05, 0.05, 0.05, 0.05] ./ t,
    ΔX_min = [0.005, 0.005, 0.005, 0.005] ./ t,
    ΔX_max = [3.0, 3.0, 1.0, 1.0] ./ t,
    ΔU_initial = [0.3, 0.1],
    ΔU_min = [0.03, 0.01],
    ΔU_max = [2.0, 0.5],
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
    # Cap funnels inside the domain box by construction — without it the
    # :maximin funnel at the arc→reverse junction balloons past the domain and
    # the state-domain gate rejects a posteriori (measured, k = 28/54). The
    # walls stay a-posteriori: the cap is per-axis slabs on the included box.
    domain_cap = true,
    check_state_domain = true,
)

println("— generate ⇄ certify loop (retry ladder, up to 3 rounds) —")
bw = EB.BackwardCertifier(zprovider, sdp, back_opts)
driver = AB.TrajectoryCertificationOptimizer.Optimizer(mppi, bw)
MOI.set(driver, MOI.RawOptimizerAttribute("concrete_problem"), gen_problem)
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

lifted = lift(traj)
# Nominal-in-domain check: wherever the nominal leaves the cover domain, the
# per-step domain cap silently disengages and the gate rejects a posteriori.
let xs = collect(ST.states(lifted)), out = findall(x -> collect(x) ∉ cover.system.X, xs)
    if !isempty(out)
        println("  WARNING: nominal outside the cover domain at state(s) $out:")
        for i in out
            println("    x[$i] = $(round.(collect(xs[i]); digits = 3))")
        end
    end
end
let xe = collect(ST.states(lifted))[end]
    r_box = 0.9 .* [0.5, 0.5, deg2rad(5.0), deg2rad(5.0)]
    margin = sum(((collect(xe) .- collect(x_goal)) ./ r_box) .^ 2)
    println(
        "  endpoint = $(round.(collect(xe); digits = 3)); box-centered margin = ",
        "$(round(margin; digits = 3)) (engages iff ≤ 0.64)",
    )
end

# ------------------------------------------------------------
# 4b) Capture-phase fallback (the double-pendulum pattern): when the chain dies
# in the benign forward arc, re-certify the certified suffix STANDALONE as a
# complete chain — the guarantee then reads "from anywhere in the entry funnel,
# the controller parks the rig", which for this benchmark covers the entire
# jackknife-unstable reverse phase.
# ------------------------------------------------------------

if !loop_success && bres.failed_k !== nothing && !isempty(bres.lmi_data.ellipsoids)
    ksuf = bres.failed_k + 1
    xs_l = collect(ST.states(lifted))
    us_l = collect(ST.inputs(lifted))
    suffix = ST.Trajectory(xs_l[ksuf:end]; inputs = us_l[ksuf:end])
    println(
        "— capture-phase certificate: re-certifying the suffix (states $ksuf..$(length(xs_l))) —",
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
# 5) Closed-loop falsification: sample a meaningful-volume funnel ellipsoid,
# replay the certified feedbacks on the plant map, count target hits.
# ------------------------------------------------------------

if !isempty(ells) && length(ells) > 1
    i_val = min(something(findfirst(v -> v >= 1e-6, vols), 1), length(ells) - 1)
    E_val = ells[i_val]
    f_z(z, u) = SVector{4}(f(SVector{4}(t .* z), u) ./ t)
    zT = zbox(T_plus)
    n_ok = 0
    n_samples = 50
    rng_val = Random.MersenneTwister(7)
    for zs0 in UT.samples(E_val, n_samples; rng = rng_val)
        z = SVector{4}(zs0)
        for κ in bres.lmi_data.kappas[i_val:end]
            u = Matrix(κ.A) * collect(z) .+ collect(κ.c)
            u = SVector(clamp(u[1], -5.0, 5.0), clamp(u[2], -pi / 4, pi / 4))
            z = f_z(z, u)
        end
        global n_ok += (z ∈ zT)
    end
    println(
        "  closed-loop from the $(round((length(ells) - i_val) * Δt; digits = 2))-s-out ",
        "region (V = $(round(vols[i_val]; sigdigits = 3))): ",
        "$(n_ok)/$(n_samples) samples reach the target",
    )
end

# ------------------------------------------------------------
# 6) Plot: parking lot with walls, trajectory, funnel shadows + bay zoom
# ------------------------------------------------------------

fig = plot(;
    xlabel = "x  [m]",
    ylabel = "y  [m]",
    title = "Certified reverse parking (tractor-trailer)",
    legend = :topleft,
    size = (900, 700),
    aspect_ratio = :equal,
)

box2(H) = LazySets.Hyperrectangle(;
    low = SVector(LazySets.low(H, 1), LazySets.low(H, 2)),
    high = SVector(LazySets.high(H, 1), LazySets.high(H, 2)),
)
AV.plot_xy_obstacles!(fig, AV.parking_obstacles(); alpha = 0.6, color = :dimgray)
plot!(fig, box2(problem.initial_set); color = :gray, alpha = 0.5, label = "initial set")
plot!(fig, box2(problem.target_set); color = :green, alpha = 0.35, label = "target (bay)")

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

xs_plot = collect(ST.states(lifted))
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
    title = "zoom: backing into the bay",
    legend = false,
    aspect_ratio = :equal,
    xlims = (6.0, 10.3),
    ylims = (4.3, 6.7),
)
AV.plot_xy_obstacles!(figz, AV.parking_obstacles(); alpha = 0.6, color = :dimgray)
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

final = plot(fig, figz; layout = (1, 2), size = (1500, 680))
plot_path = joinpath(@__DIR__, "demo_vehicle_parking.png")
savefig(final, plot_path)
println("— plot saved: $plot_path —")
