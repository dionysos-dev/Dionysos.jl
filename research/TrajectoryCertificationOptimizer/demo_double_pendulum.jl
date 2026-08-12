# Double-pendulum swing-up: MPPI generation + ellipsoidal funnel certification.
#
# Task (problems/Pendulum/double_pendulum.jl, swing_up sets): swing the FREE second
# link upright (θ₂ ∈ π ± 50°) with one torque on the shoulder (u enters θ̈₁ only —
# underactuated), θ₁ hard-limited to ±π/2, ω ∈ ±5, from both links hanging (±3°).
# The benchmark's deadband input set (set_minus ±5.5 \ ±0.5) is replaced by the
# plain ±5.5 box: the certificate must be sound against a realizable actuator
# (same decision as the simple-pendulum demo).
#
# The certified plant is the RK2 (midpoint) map at Δt = 0.025 — generation and
# certification run on the IDENTICAL map (no model mismatch at all). RK4 is not
# an option here: its four symbolic self-compositions of the rational-trig
# dynamics (1/(m₁+m₂sin²Δθ)) make the provider's Hessian expressions overflow
# the compiler stack (measured; Euler×4 substeps dies the same way — composition
# depth, not stage structure). RK2's one composition compiles in ~20 s and its
# normalized Lipschitz bounds are within 30% of Euler's.
#
# Generation (raw CEM, no abstraction seed): the hanging state is a CEM trap —
# random torque bursts exit θ₁ ∈ ±π/2 and pay the domain penalty, so doing
# nothing wins. Two additions unlock it: a dense height-shaping stage cost on
# link 2 (1 + cos θ₂ — periodic, seam-blind) and a resonant-forcing seed scan
# u = A·sin(ωt + φ) (θ₁ responds ~A/ω², staying in the tube; link 2 is driven
# near its hanging frequency √(g/l) ≈ 3.1). CEM then finds a direct ~6 s
# swing-up ending deep in the target (box-centered terminal seed engages).
#
# Certification: globally normalized coordinates, :maximin, :ball remainder,
# adaptive boxes, and the NEW domain_cap — per-step SOC slabs confining every
# funnel to the state domain AND to its linearization box by construction
# (without it, the size-maximizing SDP grows box-sized funnels that overflow
# θ₁/ω a posteriori — measured k=89 domain-gate failure — or drives the box
# search into Hessian-bound blowup). The lifted θ₂ is capped to the swing
# corridor ±1 rad: not a physical claim (θ₂ is periodic — only speeds carry
# hard bounds) but SDP conditioning + unambiguous wrap-back of the funnel
# chart (funnel θ₂-width must stay below half a period).
#
# MEASURED WALL (four-lever ablation, all at the same physical state — second
# link ~80° below horizontal on the upswing, |ω₂| ≈ 2.5–3):
#   headroom 0.5 → 1.0    : runway 27 steps (1.35 s) — unchanged
#   soft pace cap |ω₂|≲3.5 : same wall, same runway
#   Δt 0.05 → 0.025       : 27 → 57 steps = 1.35 → 1.43 s — Δt-INVARIANT
#   free-angle domain      : θ₁ freed to the full circle (only ω bounded, both
#     angles periodic/lifted) — the optimal swing never leaves ±π/2 anyway
#     (lifted θ₁ corridor [-0.81, 1.26]); 27 → 30 steps, same wall, same Vmax.
# The funnel bleed through the ballistic ascent is per-unit-time intrinsic:
# one (saturated) input defending 4 states through the most deviation-
# amplifying stretch. Certified volumes decay 7.76 → 1e-10 over ~1.4 s.
#
# SHIPPED RESULT — the capture-phase certificate (user-approved fallback): the
# certified suffix re-certified standalone as a complete chain — a ~1.4 s funnel
# corridor into the upright target, terminal ⊆ target, domain gates green. The
# corridor is reported as a DEPTH PROFILE (guaranteed recovery region vs time
# before target: ~[0.073, 0.073, 0.36, 0.36] physical semi-axes at ~0.1 s out,
# needle at the 1.4 s cliff), and the FunnelController is validated closed-loop
# on samples from the meaningful-region depth.
#
# REMAINDER-MODEL LADDER (measured on the up-up trajectory, transfers here):
# :ball and :john_ball are equivalent (both ball covers are corner-tight; the
# per-axis radii buy nothing), :vertices is strictly tighter (exact joint
# corner support). Regime dependence: at Δt = 0.05 (up-up) :vertices moves the
# WALL (+33% runway); at Δt = 0.025 (here) the per-step box is small enough
# that the wall is physics-dominated — runway unchanged, but the guaranteed
# regions DOUBLE at every depth. 16 blocks/SDP is affordable at n = 4.
#
# vs PR #569 (Florentin, run_double_pendulum_mppi) — NOT the same task: his
# MPPI experiment ran `benchmark_up_convex`, the UP-UP HANDSTAND (both angles
# in π ± 25°, torque ±8.5, ω domain ±6, both angles periodic), where his chain
# died at k=34/66 via logdet pancake collapse (radii ~3e-8) after fixed-box
# steps whose own audit shows required radii exceeding the box 14.8× near the
# terminal (consistency not enforced ⇒ those steps are informal). This demo
# ships the swing_up objective (second link up, shoulder in ±π/2); the up-up
# head-to-head with our stack is a separate experiment. His plant also
# differs: an operator-precedence slip leaves the Coriolis term undivided by
# (l·α) in BOTH his dynamic and his symbolic model (internally consistent,
# but not the textbook double pendulum this demo certifies).

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const AB = DI.Optim.Abstraction
const EB = AB.EllipsoidalTrajectoryCertifier

import LazySets
import LinearAlgebra as LA
import MathematicalSystems as MS
using StaticArrays
using Random
using Symbolics
import Clarabel
using JuMP: optimizer_with_attributes
using Plots

const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")
include(joinpath(PROBLEMS, "Pendulum", "double_pendulum.jl"))
const DP = DoublePendulum

params = DP.Params()

_X_ = LazySets.Hyperrectangle(;
    low = SVector(-π / 2, -π, -5.0, -5.0),
    high = SVector(π / 2, π, 5.0, 5.0),
)
_U_ = LazySets.Hyperrectangle(; low = SVector(-5.5), high = SVector(5.5))
_I_ = LazySets.Hyperrectangle(;
    low = SVector(-3π / 180, -3π / 180, -0.5, -0.5),
    high = SVector(3π / 180, 3π / 180, 0.5, 0.5),
)
_T_ = LazySets.Hyperrectangle(;
    low = SVector(-π / 2, π - 50π / 180, -4.5, -4.5),
    high = SVector(π / 2, π + 50π / 180, 4.5, 4.5),
)

Δt = 0.025
f_ode = DP.dynamic(params)
# The RK2 plant map — the ONE map generated on, certified, and simulated.
f = (x, u) -> x .+ Δt .* f_ode(x .+ (Δt / 2) .* f_ode(x, u), u)
dsys = MS.ConstrainedBlackBoxControlDiscreteSystem(f, 4, 1, _X_, _U_)

# θ₂ is the periodic dimension (θ₁ is hard-limited by X, not a seam).
wrap = UT.get_periodic_wrapper(SVector(2), SVector(2π); start = SVector(-π))
T_split = UT.set_in_period(_T_, SVector(2), SVector(2π), SVector(-π))
discrete_problem =
    PR.OptimalControlProblem(dsys, _I_, T_split, nothing, nothing, PR.Infinity(), nothing)

# ------------------------------------------------------------
# 1) Generation: resonant seed scan + height-shaped CEM
# ------------------------------------------------------------

println("— MPPI (CEM) —")
x0 = SVector(0.0, 0.0, 0.0, 0.0)
nstep = 240
u_clamp = 4.5    # reserve 1.0 vs the plant's ±5.5 (certification headroom)

function rollout_inputs(us)
    xs = [x0]
    for u in us
        push!(xs, f(xs[end], u))
    end
    return ST.Trajectory(xs; inputs = collect(us))
end

# Wrap-aware terminal pull in the target's inscribed-ellipsoid metric.
rT = 0.9 .* [π / 2, 50π / 180, 4.5, 4.5]
struct DPTerminalPull <: AB.AbstractCostTerm
    w::Float64
end
function AB.cost_final(term::DPTerminalPull, acc, xT)
    xw = wrap(xT)
    dθ2 = rem(xw[2] - π, 2π, RoundNearest)
    return acc +
           term.w *
           ((xw[1] / rT[1])^2 + (dθ2 / rT[2])^2 + (xw[3] / rT[3])^2 + (xw[4] / rT[4])^2)
end

# Dense swing-up shaping: height of link 2 (1 + cos θ₂ ∈ [0,2], zero when up).
struct HeightShapingCost <: AB.AbstractCostTerm
    w::Float64
end
AB.cost_step(term::HeightShapingCost, acc, x, u, k) = acc + term.w * (1 + cos(x[2]))

# Soft pace cap on the second link (documented ablation: does not move the
# certification wall, kept for the gentler — direct, no wind-through — swing).
struct SwingPaceCost <: AB.AbstractCostTerm
    w::Float64
    ω_soft::Float64
end
AB.cost_step(term::SwingPaceCost, acc, x, u, k) =
    acc + term.w * max(0.0, abs(x[4]) - term.ω_soft)^2

score = AB.CompositeCost(
    AB.ReachObjectiveCost(T_split; wrap = wrap),
    DPTerminalPull(500.0),
    HeightShapingCost(5.0),
    SwingPaceCost(3.0, 3.5),
    AB.InputEffortCost(0.001),
    AB.InputSmoothnessCost(; w_du = 0.02, w_ddu = 0.005),
    AB.DomainPenaltyCost(_X_; wrap = wrap),
)

# Resonant-forcing seed scan: u = A·sin(ωt + φ), scored with the CEM cost.
import Dionysos.Optim.Abstraction: rollout_cost
function scan_seeds()
    best_seed, best_cost, best_cfg = nothing, Inf, nothing
    for A in (2.0, 3.0, 4.0, 5.0), ω in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0), φ in (0.0, π / 2, π)
        us = [
            SVector(clamp(A * sin(ω * k * Δt + φ), -u_clamp, u_clamp)) for
            k in 0:(nstep - 1)
        ]
        c, _ = rollout_cost(f, x0, us, wrap, score, _X_, false, 1e6)
        if c < best_cost
            best_cost, best_cfg = c, (A, ω, φ)
            best_seed = rollout_inputs(us)
        end
    end
    return best_seed, best_cost, best_cfg
end
seed_traj, seed_cost, seed_cfg = scan_seeds()
println("  seed scan: cost = $(round(seed_cost; digits = 1)), (A, ω, φ) = $seed_cfg")

mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = seed_traj,
    nstep = nstep,
    nsamples = 3000,
    niter = 120,
    anneal = 0.99,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.8)),
    project_input = u -> SVector(clamp(u[1], -u_clamp, u_clamp)),
    cost = score,
    wrap_state = (p, x) -> wrap(x),
    update_rule = :cem,
    elite_frac = 0.05,
    antithetic = true,
    # Keep optimizing past the first success: the pull terms center the endpoint,
    # which is what lets the box-centered terminal seed engage.
    stop_on_success = false,
)
AB.set_problem!(mppi, discrete_problem)
gen_time = @elapsed AB.generate!(mppi)
traj = AB.get_trajectory(mppi)
println(
    "  $(round(gen_time; digits = 1)) s, success = $(AB.get_success(mppi)), ",
    "steps = $(length(ST.inputs(traj)))",
)
AB.get_success(mppi) || error("generation failed — retune before certifying")

# Terminal-seed diagnostics (endpoint depth in the target).
let xe = wrap(collect(ST.states(traj))[end])
    c2 = rem(xe[2] - π, 2π, RoundNearest)
    margin = (xe[1] / rT[1])^2 + (c2 / rT[2])^2 + (xe[3] / rT[3])^2 + (xe[4] / rT[4])^2
    println(
        "  endpoint (wrapped) = $(round.(collect(xe); digits = 3)), ",
        "box-centered margin = $(round(margin; digits = 3)) (engages iff ≤ 0.64)",
    )
end

# ------------------------------------------------------------
# 2) x-frame lift: θ₂ periodic unwrap, shifted so the endpoint lands near π
# ------------------------------------------------------------

lift = function (tr)
    lifted = ST.unwrap_trajectory(tr, (2,), (2π,))
    θN = collect(ST.states(lifted))[end][2]
    shift = 2π * round((θN - π) / (2π))
    shift == 0.0 && return lifted
    xs = [SVector(x[1], x[2] - shift, x[3], x[4]) for x in ST.states(lifted)]
    return ST.Trajectory(xs; inputs = collect(ST.inputs(lifted)))
end
lifted_traj = lift(traj)
θ2s = [x[2] for x in ST.states(lifted_traj)]
println(
    "  lifted θ₂ corridor: [$(round(minimum(θ2s); digits = 2)), ",
    "$(round(maximum(θ2s); digits = 2))]",
)

# ------------------------------------------------------------
# 3) Symbolic provider on the SAME RK2 map, then normalized coordinates
# ------------------------------------------------------------

Symbolics.@variables xv[1:4] uv[1:1] wv[1:4]
xsym = collect(xv)
usym = collect(uv)
wsym = collect(wv)
f_cont(xloc, uloc) = collect(DP.dynamic(params)(xloc, uloc))
rk2 = xsym .+ Δt .* f_cont(xsym .+ (Δt / 2) .* f_cont(xsym, usym), usym)
fsymbolic = [rk2[i] + wsym[i] for i in 1:4]
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
    ST.format_input_set(_U_),                # the plant's true (convex) input set
    ST.format_noise_set(Wset),
)

t = [0.15, 0.15, 0.75, 0.75]                 # characteristic scales (rad, rad, rad/s, rad/s)
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{4}(LazySets.low(H) ./ t),
    high = SVector{4}(LazySets.high(H) ./ t),
)
# Certified region: real bounds on θ₁ (spec joint limit) and ω (domain), the
# lifted θ₂ bounded to the swing corridor (conditioning + wrap-back chart width —
# NOT a physical wall; θ₂ is periodic).
X_gate = LazySets.Hyperrectangle(;
    low = SVector(-π / 2, minimum(θ2s) - 1.0, -5.0, -5.0),
    high = SVector(π / 2, maximum(θ2s) + 1.0, 5.0, 5.0),
)
zsys = MS.ConstrainedBlackBoxControlContinuousSystem(f_ode, 4, 1, zbox(X_gate), _U_)
zproblem = PR.OptimalControlProblem(
    zsys,
    zbox(_I_),
    zbox(_T_),
    nothing,
    nothing,
    PR.Infinity(),
    nothing,
)
ztraj(tr) = ST.Trajectory(
    [SVector{4}(collect(x) ./ t) for x in ST.states(tr)];
    inputs = collect(ST.inputs(tr)),
)
ztr = ztraj(lifted_traj)

# ------------------------------------------------------------
# 4) Backward certification: domain-capped, :maximin, :ball, box-scale ladder
# ------------------------------------------------------------

sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

# With domain_cap the funnel is confined to each candidate box, so the scale
# ladder IS the funnel-size search.
adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    true,
    [0.05, 0.05, 0.05, 0.05] ./ t,
    [0.005, 0.005, 0.005, 0.005] ./ t,
    [1.5, 1.5, 4.0, 4.0] ./ t,
    [0.2],
    [0.02],
    [1.0],
    1.5,
    1.05,
    30,
    1e-8,
    false,
    [1.0, 2.0, 4.0, 8.0, 16.0, 32.0],
    :max_volume,
    false,
)
back_opts = EB.ChainOptions(;
    maxδx = 8.0,
    maxδu = 1.0,
    λ = 0.001,
    terminal_shape = nothing,
    terminal_shrink = 0.9,
    state_scaling = nothing,                 # exact: dynamics already normalized
    linearization_δx = [0.05, 0.05, 0.05, 0.05] ./ t,
    linearization_δu = [0.2],
    adaptive_boxes = adaptive_opts,
    objective = :maximin,
    # Measured remainder ladder on the up-up trajectory (same chain machinery):
    # :ball wall k=32, :john_ball k=32 (per-axis radii buy nothing — both balls
    # are corner-tight; the win is the exact joint corner support), :vertices
    # k=24 (+33% runway). 2⁴ = 16 blocks/SDP is affordable at this chain length.
    remainder_model = :vertices,
    domain_cap = true,
    check_state_domain = true,
)

println("— full-chain certification attempt —")
bw = EB.BackwardCertifier(zprovider, sdp, back_opts)
EB.set_problem!(bw, zproblem)
EB.set_trajectory!(bw, ztr)
cert_time = @elapsed EB.certify!(bw)
bres = EB.get_result(bw)
println(
    "  $(round(cert_time; digits = 1)) s, certified = $(bres.success), ",
    "failed_k = $(bres.failed_k), terminal ⊆ target: $(bres.terminal_contained), ",
    "domain gate: $(bres.state_domain_checked)",
)

# ------------------------------------------------------------
# 5) The capture-phase certificate: re-certify the certified suffix standalone
# ------------------------------------------------------------

zxs = collect(ST.states(ztr))
zus = collect(ST.inputs(ztr))
if bres.success
    capture_bw = bw
    capture_start = 1
else
    capture_start = bres.failed_k + 1
    println("— capture chain: re-certifying the suffix from state $capture_start —")
    capture_ztr = ST.Trajectory(zxs[capture_start:end]; inputs = zus[capture_start:end])
    capture_bw = EB.BackwardCertifier(zprovider, sdp, back_opts)
    EB.set_problem!(capture_bw, zproblem)
    EB.set_trajectory!(capture_bw, capture_ztr)
    capture_time = @elapsed EB.certify!(capture_bw)
    cres = EB.get_result(capture_bw)
    println(
        "  $(round(capture_time; digits = 1)) s, certified = $(cres.success), ",
        "steps = $(length(zus) - capture_start + 1) ",
        "($(round((length(zus) - capture_start + 1) * Δt; digits = 2)) s of trajectory), ",
        "terminal ⊆ target: $(cres.terminal_contained)",
    )
    cres.success || error("the capture chain must certify — it just did as a suffix")
end
cres = EB.get_result(capture_bw)

# The corridor is funnel-shaped: fat at the target, needle-thin at the cliff
# (measured: volume 14.9 → 1e-12 over 57 steps). One chain carries every depth —
# report the guaranteed recovery region vs time-before-target instead of a single
# entry number (a 1e-3 volume floor would keep only 5 steps, semi-axes
# [0.06, 0.06, 0.3, 0.3] — the initial-box scale; the full 57-step corridor ends
# in a needle. Both facts below).
ells = cres.lmi_data.ellipsoids
vols = [π^2 / 2 * sqrt(LA.det(D * Matrix(LazySets.shape_matrix(E)) * D)) for E in ells]
phys_radii(E) = t .* sqrt.(LA.eigvals(LA.Symmetric(Matrix(LazySets.shape_matrix(E)))))
println("  capture depth profile (guaranteed recovery region vs time before target):")
for d in unique(clamp.([5, 10, 20, 40, length(ells) - 1], 1, length(ells) - 1))
    E = ells[end - d]
    r = sort(phys_radii(E))
    println(
        "    $(lpad(round(d * Δt; digits = 2), 5)) s out: V = ",
        "$(round(vols[end - d]; sigdigits = 3)), semi-axes $(round.(r; sigdigits = 2))",
    )
end

zctrl = EB.get_controller(capture_bw)
println("— certified controller: FunnelController with $(length(zctrl.kappas)) steps —")

# ------------------------------------------------------------
# 6) Closed-loop validation from the MEANINGFUL-region depth: sample the last
# funnel with volume ≥ 1e-3 (a real region, not the cliff needle) and run the
# remaining controller steps.
# ------------------------------------------------------------

i_val = min(something(findfirst(v -> v >= 1e-3, vols), 1), length(ells) - 1)
E_val = ells[i_val]
f_z(z, u) = SVector{4}(f(SVector{4}(t .* z), u) ./ t)
zT = zbox(_T_)
n_ok = 0
n_samples = 50
rng_val = Random.MersenneTwister(7)
for zs0 in UT.samples(E_val, n_samples; rng = rng_val)
    z = SVector{4}(zs0)
    for κ in cres.lmi_data.kappas[i_val:end]
        u = SVector{1}(clamp.(Matrix(κ.A) * collect(z) .+ collect(κ.c), -5.5, 5.5))
        z = f_z(z, u)
    end
    global n_ok += (z ∈ zT)
end
println(
    "  closed-loop from the $(round((length(ells) - i_val) * Δt; digits = 2))-s-out ",
    "region (V = $(round(vols[i_val]; sigdigits = 3))): ",
    "$(n_ok)/$(n_samples) samples reach the target",
)

# ------------------------------------------------------------
# 7) Plots: swing plane (θ₂, ω₂) + shoulder plane (θ₁, ω₁), with funnel shadows
# ------------------------------------------------------------

denorm_shadow(E, i, j) = begin
    Q = D * Matrix(LazySets.shape_matrix(E)) * D
    c = t .* collect(LazySets.center(E))
    LazySets.Ellipsoid(c[[i, j]], Matrix(LA.Symmetric(Q[[i, j], [i, j]])))
end

xs_lift = collect(ST.states(lifted_traj))
cap_xs = xs_lift[capture_start:end]

fig1 = plot(;
    xlabel = "θ₂  [rad] (lifted)",
    ylabel = "ω₂  [rad/s]",
    title = "Double pendulum swing-up: certified capture corridor",
    legend = :topleft,
    size = (900, 650),
)
Tbox2 = LazySets.Hyperrectangle(;
    low = SVector(π - 50π / 180, -4.5),
    high = SVector(π + 50π / 180, 4.5),
)
plot!(fig1, Tbox2; color = :green, alpha = 0.3, label = "target (θ₂, ω₂)")
for (i, E) in enumerate(ells)
    plot!(
        fig1,
        denorm_shadow(E, 2, 4);
        color = :steelblue,
        alpha = 0.25,
        linewidth = 1.0,
        linecolor = :steelblue,
        label = i == 1 ? "funnel (θ₂ω₂ shadow)" : "",
    )
end
plot!(
    fig1,
    [x[2] for x in xs_lift],
    [x[4] for x in xs_lift];
    color = :gray,
    linewidth = 1.5,
    label = "nominal (full swing)",
)
plot!(
    fig1,
    [x[2] for x in cap_xs],
    [x[4] for x in cap_xs];
    color = :black,
    linewidth = 2.5,
    label = "certified capture phase",
)

fig2 = plot(;
    xlabel = "θ₁  [rad]",
    ylabel = "ω₁  [rad/s]",
    title = "shoulder plane",
    legend = :topleft,
)
plot!(
    fig2,
    LazySets.Hyperrectangle(; low = SVector(-π / 2, -4.5), high = SVector(π / 2, 4.5));
    color = :green,
    alpha = 0.3,
    label = "target (θ₁, ω₁)",
)
for E in ells
    plot!(
        fig2,
        denorm_shadow(E, 1, 3);
        color = :steelblue,
        alpha = 0.25,
        linewidth = 1.0,
        linecolor = :steelblue,
        label = "",
    )
end
plot!(
    fig2,
    [x[1] for x in xs_lift],
    [x[3] for x in xs_lift];
    color = :gray,
    linewidth = 1.5,
    label = "nominal",
)
plot!(
    fig2,
    [x[1] for x in cap_xs],
    [x[3] for x in cap_xs];
    color = :black,
    linewidth = 2.5,
    label = "capture phase",
)

final_fig = plot(fig1, fig2; layout = (1, 2), size = (1500, 650))
plot_path = joinpath(@__DIR__, "demo_double_pendulum.png")
savefig(final_fig, plot_path)
println("— plot saved: $plot_path —")
