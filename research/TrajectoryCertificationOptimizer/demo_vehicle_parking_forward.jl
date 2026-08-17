# Certified forward (nose-in) parking — the same "Permis" lot as
# demo_vehicle_parking.jl, entered head-first.
#
# Same 4-D tractor-trailer, same compact rig (L1 = 1, L2 = 1, Lc = 0.5), same
# two walls forming the dead-end bay (corridor y ∈ [4.7, 6.0] for x ∈ [4, 10]).
# The vehicle starts near the origin heading east and drives INTO the bay
# nose-first through a forward S (θ stays near 0 — no U-turn, no periodic
# cover, no unwrap). The forward direction is the ϕ-STABLE one: the trailer
# follows instead of jackknifing, which is exactly why this maneuver is the
# easy sibling of the reverse bay demo — same corridor, opposite dynamics
# regime.
#
#     julia --project=test research/TrajectoryCertificationOptimizer/demo_vehicle_parking_forward.jl
#
# Pipeline: scripted forward-S seed scan (arc + counter-arc + heading-feedback
# cruise, domain-clean scripts only) → CEM-MPPI → backward :maximin
# certification in normalized coordinates with the state-domain gate ACTIVE
# against the walls → capture-phase fallback → closed-loop falsification.

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
problem = AV.parking_problem(; params = params, heading = :nose_in)

Δt = 0.2
base = PR.discretize_problem(problem, Δt; num_substeps = 4)
f = MS.mapping(base.system)

# ------------------------------------------------------------
# 1) Seed: forward-S script scan (arc up, counter-arc east, cruise into the bay)
# ------------------------------------------------------------

x0 = SVector(0.0, 0.0, 0.0, 0.0)
x_goal = SVector(8.35, 5.5, 0.0, 0.0)          # center of the nose-in target

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
    # Cruise east into the bay under a heading feedback (forward driving is
    # ϕ-stable; only θ needs steering).
    for _ in 1:n3
        x = xs[end]
        step!(SVector(2.0, clamp(-1.5 * x[3], -0.45, 0.45)))
    end
    return xs, us
end

score(x) = (x[1] - x_goal[1])^2 + (x[2] - x_goal[2])^2 + 4 * x[3]^2 + 4 * x[4]^2

best = (Inf, 0.4, 8, 0.4, 8, 10)
for δ1 in 0.2:0.05:0.6, n1 in 4:2:16, δ2 in 0.2:0.05:0.6, n2 in 4:2:16, n3 in 0:4:20
    xs, _ = rollout_script(δ1, n1, δ2, n2, n3)
    all(x -> collect(x) ∈ problem.system.X, xs) || continue
    s = score(xs[end])
    s < best[1] && (global best = (s, δ1, n1, δ2, n2, n3))
end
isfinite(best[1]) || error("seed scan: no domain-clean script found")
_, δ1, n1, δ2, n2, n3 = best
println(
    "— seed scan: arc δ = $δ1 × $n1, counter-arc δ = $δ2 × $n2, cruise × $n3, " *
    "endpoint score = $(round(best[1]; digits = 3)) —",
)

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
u_plan = SVector(2.5, 0.45)                    # reserve vs the plant's [±5, ±π/4]

cost = AB.CompositeCost(
    AB.ReachObjectiveCost(problem.target_set),
    AB.TerminalPullCost(collect(x_goal), [0.3, 0.45, deg2rad(4.5), deg2rad(4.5)]; w = 1e4),
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
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.15, 0.06)),
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

t = [2.0, 2.0, 0.35, 0.35]                     # characteristic scales (m, m, rad, rad)
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{4}(LazySets.low(H) ./ t),
    high = SVector{4}(LazySets.high(H) ./ t),
)

lot_box = LazySets.Hyperrectangle(;
    low = SVector(-1.0, -1.0, -pi, -deg2rad(50.0)),
    high = SVector(10.0, 9.0, pi, deg2rad(50.0)),
)
zX = UT.set_minus(
    zbox(lot_box),
    UT.set_union([
        zbox(AV.extrude_xy_obstacle_to_4d(ob, lot_box)) for ob in AV.parking_obstacles()
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
    r_box = 0.9 .* [0.35, 0.5, deg2rad(5.0), deg2rad(5.0)]
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
            SVector(clamp(u[1], -5.0, 5.0), clamp(u[2], -pi / 4, pi / 4)),
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
# 6) Plot: the lot with walls, trajectory, funnel shadows + bay zoom
# ------------------------------------------------------------

fig = plot(;
    xlabel = "x  [m]",
    ylabel = "y  [m]",
    title = "Certified forward parking (nose-in)",
    legend = :topleft,
    size = (900, 700),
    aspect_ratio = :equal,
)

box2(H) = LazySets.Hyperrectangle(;
    low = SVector(LazySets.low(H, 1), LazySets.low(H, 2)),
    high = SVector(LazySets.high(H, 1), LazySets.high(H, 2)),
)
AV.plot_xy_obstacles!(fig, AV.parking_obstacles())
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
    title = "zoom: driving into the bay",
    legend = false,
    aspect_ratio = :equal,
    xlims = (3.0, 10.2),
    ylims = (4.0, 7.0),
)
AV.plot_xy_obstacles!(figz, AV.parking_obstacles())
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
plot_path = joinpath(@__DIR__, "demo_vehicle_parking_forward.png")
savefig(final, plot_path)
println("— plot saved: $plot_path —")

# ------------------------------------------------------------
# 7) Dashboard animation: rig in the parking lot + xy and input panels
# ------------------------------------------------------------

av_plot = AV.system_plot!(;
    params = params,
    obstacles2d = AV.parking_obstacles(),
    xlims = (-4.0, 10.5),
    ylims = (-2.0, 9.5),
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
    filename = joinpath(@__DIR__, "demo_vehicle_parking_forward_dashboard.gif"),
    title = "Certified forward parking",
)
println("— dashboard saved: $dash_path —")
