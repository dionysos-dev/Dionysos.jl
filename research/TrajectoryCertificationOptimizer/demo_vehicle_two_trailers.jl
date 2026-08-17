# Certified two-trailer corridor parking — the 5-D model from PR #569, ported
# to the modern API and given the benchmark its placeholder factory never had.
#
# 5-D tractor + TWO trailers (x, y, θ, ϕ₁, ϕ₂), rig L1 = 2.2, L2 = 2.8,
# L3 = 2.2, Lc = Lc2 = 0.9 — ~10 m nose-to-tail. Two walls leave a dead-end
# corridor (y ∈ [4.5, 7.5] for x ∈ [12, 30]) opening west; the rig starts in
# the open area heading east and drives nose-first into the corridor, parking
# deep inside with both hitch angles settled. Forward is the ϕ-STABLE
# direction — both trailers follow instead of jackknifing — which is what
# makes the 5-D chain certifiable end-to-end.
#
#     julia --project=test research/TrajectoryCertificationOptimizer/demo_vehicle_two_trailers.jl
#
# Pipeline: scripted forward-S seed scan (arc + counter-arc + heading-feedback
# cruise, domain-clean scripts only) → CEM-MPPI → backward :maximin
# certification in normalized coordinates with the :ball remainder (n = 5:
# the measured ladder rule — :vertices' 2ⁿ corner blocks stop paying past
# n = 4) and the state-domain gate ACTIVE against the walls → capture-phase
# fallback → closed-loop falsification + inflation stress.
#
# Status (2026-08-17, Clarabel, fixed rng): FULL-CHAIN certificate, round 1 —
# the entire 55-step maneuver certifies end-to-end (failed_k = nothing,
# terminal ⊆ target, wall-aware domain gates green), endpoint
# (24.00, 6.00, 0, 0, 0) = the target center exactly (box-centered terminal
# engaged, Vmax ≈ 1e-3), the domain walls carrying a 0.6 m body-clearance
# margin so the drawn ~1 m rig stays visually clear of the TRUE walls.
# Inflation stress: 100% empirical success out to α = 3 — the forward
# ϕ-stable direction gives the widest recovery basin of the vehicle family.
# COST: the 5-D LMI chain is expensive (~80 min for the chain: bigger PSD
# blocks × adaptive box iterations) — this demo is a long run, not an
# iteration loop. The 5-D model itself (PR #569's one genuinely new system)
# was never executed on that branch; this is its first certified run.

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
include(joinpath(PROBLEMS, "ArticulatedVehicle", "articulated_vehicle_2trailers.jl"))
const AV2 = ArticulatedVehicle2Trailers

params = AV2.Params()
# Lateral body-clearance margin: the domain walls are inflated by ~half the
# rig's width so the certified POINT trajectory keeps the drawn bodies clear
# of the TRUE walls (which the plots and dashboard show).
obstacle_margin = 0.6
problem = AV2.problem(;
    params = params,
    obstacles2d = AV2.corridor_obstacles(; margin = obstacle_margin),
)

Δt = 0.25
base = PR.discretize_problem(problem, Δt; num_substeps = 4)
f = MS.mapping(base.system)

# ------------------------------------------------------------
# 1) Seed: forward-S script scan (arc up, counter-arc east, cruise inside)
# ------------------------------------------------------------

x0 = SVector(2.5, 2.0, 0.0, 0.0, 0.0)
x_goal = SVector(24.0, 6.0, 0.0, 0.0, 0.0)     # center of the corridor target

function rollout_script(δ1, n1, δ2, n2, n3)
    xs = [x0]
    us = SVector{2, Float64}[]
    step!(u) = (push!(us, u); push!(xs, f(xs[end], u)))
    for _ in 1:n1
        step!(SVector(2.0, δ1))
    end
    for _ in 1:n2
        step!(SVector(2.0, -δ2))
    end
    # Cruise east under a heading feedback (forward driving is ϕ-stable for
    # both trailers; only θ needs steering).
    for _ in 1:n3
        x = xs[end]
        step!(SVector(2.0, clamp(-1.5 * x[3], -0.5, 0.5)))
    end
    return xs, us
end

score(x) =
    (x[1] - x_goal[1])^2 + (x[2] - x_goal[2])^2 + 4 * x[3]^2 + 4 * x[4]^2 + 4 * x[5]^2

println("— seed scan starting —")
flush(stdout)
best = (Inf, 0.4, 8, 0.4, 8, 20)
scan_time = @elapsed for δ1 in 0.15:0.05:0.5,
    n1 in 4:2:16,
    δ2 in 0.15:0.05:0.5,
    n2 in 4:2:16,
    n3 in 8:4:36

    xs, _ = rollout_script(δ1, n1, δ2, n2, n3)
    all(x -> collect(x) ∈ problem.system.X, xs) || continue
    s = score(xs[end])
    s < best[1] && (global best = (s, δ1, n1, δ2, n2, n3))
end
println("  scan: $(round(scan_time; digits = 1)) s")
flush(stdout)
isfinite(best[1]) || error("seed scan: no domain-clean script found")
_, δ1, n1, δ2, n2, n3 = best
println(
    "— seed scan: arc δ = $δ1 × $n1, counter-arc δ = $δ2 × $n2, cruise × $n3, " *
    "endpoint score = $(round(best[1]; digits = 3)) —",
)
flush(stdout)

nstep = n1 + n2 + n3 + 10                      # slack for CEM to reshape
seed_traj = begin
    xs, us = rollout_script(δ1, n1, δ2, n2, n3)
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
flush(stdout)
u_plan = SVector(2.2, 0.5)                     # reserve vs the plant's [±3, ±0.7]

cost = AB.CompositeCost(
    AB.ReachObjectiveCost(problem.target_set),
    AB.TerminalPullCost(
        collect(x_goal),
        [0.9, 0.45, deg2rad(4.5), deg2rad(4.5), deg2rad(4.5)];
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
    niter = 120,
    anneal = 0.99,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.15, 0.05)),
    project_input = u ->
        SVector(clamp(u[1], -u_plan[1], u_plan[1]), clamp(u[2], -u_plan[2], u_plan[2])),
    cost = cost,
    update_rule = :cem,
    elite_frac = 0.05,
    antithetic = true,
    stop_on_success = false,
)

# ------------------------------------------------------------
# 3) Normalized symbolic provider (5-D)
# ------------------------------------------------------------

Symbolics.@variables xv[1:5] uv[1:2] wv[1:5] T
xsym = collect(xv)
usym = collect(uv)
wsym = collect(wv)
f_cont(xloc, uloc) = collect(AV2.dynamic(params)(xloc, uloc))
f_disc = ST.runge_kutta4(f_cont, xsym, usym, T, 1)
fsymbolic = Symbolics.substitute([f_disc[i] + wsym[i] for i in 1:5], Dict(T => Δt))
Wset = LazySets.Hyperrectangle(;
    low = zero(SVector{5, Float64}),
    high = zero(SVector{5, Float64}),
)
provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    xsym,
    usym,
    wsym,
    zeros(5),
    ST.format_input_set(problem.system.U),     # the plant's true input set
    ST.format_noise_set(Wset),
)

t = [3.0, 2.0, 0.35, 0.35, 0.35]               # characteristic scales
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{5}(LazySets.low(H) ./ t),
    high = SVector{5}(LazySets.high(H) ./ t),
)

lot_box = LazySets.Hyperrectangle(;
    low = SVector(0.0, 0.0, -pi, -deg2rad(50.0), -deg2rad(50.0)),
    high = SVector(30.0, 12.0, pi, deg2rad(50.0), deg2rad(50.0)),
)
zX = UT.set_minus(
    zbox(lot_box),
    UT.set_union([
        zbox(AV2.extrude_xy_obstacle_to_5d(ob, lot_box)) for
        ob in AV2.corridor_obstacles(; margin = obstacle_margin)
    ]),
)
zsys = MS.ConstrainedBlackBoxControlContinuousSystem(
    AV2.dynamic(params),
    5,
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
# 4) Backward certification, :ball remainder, domain gate ON
# ------------------------------------------------------------

sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

adaptive_opts = EB.AdaptiveLinearizationBoxOptions(;
    ΔX_initial = [0.05, 0.05, 0.05, 0.05, 0.05] ./ t,
    ΔX_min = [0.005, 0.005, 0.005, 0.005, 0.005] ./ t,
    ΔX_max = [3.0, 3.0, 1.0, 1.0, 1.0] ./ t,
    ΔU_initial = [0.3, 0.1],
    ΔU_min = [0.03, 0.01],
    ΔU_max = [2.0, 0.5],
)
back_opts = EB.ChainOptions(;
    maxδx = 8.0,
    maxδu = 2.5,
    λ = 0.001,
    terminal_shape = nothing,
    terminal_shrink = 0.9,
    linearization_δx = [0.05, 0.05, 0.05, 0.05, 0.05] ./ t,
    linearization_δu = [0.3, 0.1],
    adaptive_boxes = adaptive_opts,
    objective = :maximin,
    remainder_model = :ball,                   # n = 5: 2ⁿ vertex blocks stop paying
    domain_cap = true,
    check_state_domain = true,
)

println("— generate ⇄ certify loop (retry ladder, up to 3 rounds) —")
flush(stdout)
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
    r_box = 0.9 .* [1.0, 0.5, deg2rad(5.0), deg2rad(5.0), deg2rad(5.0)]
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
    # 5-D ellipsoid volume: V = (8π²/15)·√det Q (the unit-5-ball volume).
    c5 = 8 * pi^2 / 15
    vols = [c5 * sqrt(LA.det(D * Matrix(LazySets.shape_matrix(E)) * D)) for E in ells]
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
    i_val = min(something(findfirst(v -> v >= 1e-8, vols), 1), length(ells) - 1)
    f_z(z, u) = SVector{5}(f(SVector{5}(t .* z), u) ./ t)
    rows = EB.inflation_stress(
        f_z,
        bres.lmi_data.kappas[i_val:end],
        bres.lmi_data.ellipsoids[i_val:end],
        zbox(problem.target_set);
        alphas = [1.0, 1.5, 2.0, 3.0],
        n_samples = 100,
        rng = Random.MersenneTwister(7),
        input_set = problem.system.U,
        project_input = u -> SVector(clamp(u[1], -3.0, 3.0), clamp(u[2], -0.7, 0.7)),
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
# 6) Plot: the lot with walls, trajectory, funnel shadows + corridor zoom
# ------------------------------------------------------------

fig = plot(;
    xlabel = "x  [m]",
    ylabel = "y  [m]",
    title = "Certified two-trailer corridor parking",
    legend = :topleft,
    size = (1000, 480),
    aspect_ratio = :equal,
    xlims = (-0.3, 30.3),
    ylims = (-0.3, 12.3),
)

box2(H) = LazySets.Hyperrectangle(;
    low = SVector(LazySets.low(H, 1), LazySets.low(H, 2)),
    high = SVector(LazySets.high(H, 1), LazySets.high(H, 2)),
)
AV2.plot_xy_obstacles!(fig, AV2.corridor_obstacles())
plot!(fig, box2(problem.initial_set); color = :gray, alpha = 0.5, label = "initial set")
plot!(
    fig,
    box2(problem.target_set);
    color = :green,
    alpha = 0.35,
    label = "target (corridor)",
)

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
    title = "zoom: entering the corridor",
    legend = false,
    aspect_ratio = :equal,
    xlims = (10.0, 28.0),
    ylims = (3.5, 8.5),
)
AV2.plot_xy_obstacles!(figz, AV2.corridor_obstacles())
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

final = plot(fig, figz; layout = (1, 2), size = (1500, 560))
plot_path = joinpath(@__DIR__, "demo_vehicle_two_trailers.png")
savefig(final, plot_path)
println("— plot saved: $plot_path —")

# ------------------------------------------------------------
# 7) Dashboard animation: three-body rig + xy and input panels
# ------------------------------------------------------------

av2_plot = AV2.system_plot!(;
    params = params,
    obstacles2d = AV2.corridor_obstacles(),
    xlims = (-0.5, 30.5),
    ylims = (-0.5, 12.5),
)
dash_plot! =
    (fg, xk, uk) -> begin
        plot!(fg, box2(problem.initial_set); color = :gray, alpha = 0.4)
        plot!(fg, box2(problem.target_set); color = :green, alpha = 0.4)
        av2_plot(fg, xk, uk)
    end
dash_path = DI.animate_trajectory_dashboard(
    dash_plot!,
    traj;
    xdims = (1, 2),
    udims = (1, 2),
    Δt = Δt,
    fps = 12,
    filename = joinpath(@__DIR__, "demo_vehicle_two_trailers_dashboard.gif"),
    title = "Certified two-trailer parking",
)
println("— dashboard saved: $dash_path —")
