# Double-pendulum UP-UP HANDSTAND: Florentin's `benchmark_up_convex` — the task
# his MPPI experiment actually ran (see his ellipsoids_rollout.gif) — with our
# full stack. Both links inverted: target θ₁, θ₂ ∈ π ± 25° (ω ∈ ±5), domain =
# full circle × ω ∈ ±6 (both angles periodic — no seam anywhere), torque ±8.5
# (plain convex box, his own addition), start hanging (±10°).
#
# Generation (raw CEM, no abstraction seed): the swing_up recipe transfers with
# three deltas — the height shaping lifts BOTH links (2 + cos θ₁ + cos θ₂), NO
# pace cost (the handstand needs the violent flip: ω rides near the ±6 domain
# wall), and a SHORT aggressive horizon (3 s; a 6 s horizon dilutes CEM and
# fails — measured). Finds a 56-step (2.8 s) double inversion, endpoint
# [3.01, 2.85, 1.45, 1.23] — both angles inside π ± 25° (his flip: 2.55 s).
#
# Certification (same machinery as demo_double_pendulum.jl: RK2 certified
# plant, normalized coordinates, :maximin, :ball, domain_cap, box-scale
# ladder): 24 enforced transitions (k=56→33), V_max ≈ 1.9 at the terminal,
# wall mid-flip at k=32 — the same per-unit-time ballistic funnel bleed as the
# swing_up (certified runway 1.2 s vs 1.4 s there). The capture-phase
# certificate is re-derived standalone below with the depth profile and
# closed-loop validation.
#
# Head-to-head on THIS task: his chain reports 32 transitions of 66, but his
# own audit shows required box radii exceeding the fixed boxes 14.8× on the
# near-terminal steps (exactly the fat, valuable ellipsoids — consistency not
# enforced ⇒ informal) before a logdet pancake collapse to radii ~3e-8; and
# his plant carries the Coriolis precedence slip. Ours: 24 transitions, every
# one consistency-enforced, collapse-proofed, domain-capped, on the textbook
# plant, with a sound fat terminal region.
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
    low = SVector(-π, -π, -6.0, -6.0),
    high = SVector(π, π, 6.0, 6.0),
)
_U_ = LazySets.Hyperrectangle(; low = SVector(-8.5), high = SVector(8.5))
_I_ = LazySets.Hyperrectangle(;
    low = SVector(-10π / 180, -10π / 180, -0.5, -0.5),
    high = SVector(10π / 180, 10π / 180, 0.5, 0.5),
)
_T_ = LazySets.Hyperrectangle(;
    low = SVector(π - 25π / 180, π - 25π / 180, -5.0, -5.0),
    high = SVector(π + 25π / 180, π + 25π / 180, 5.0, 5.0),
)

Δt = 0.05
f_ode = DP.dynamic(params)
f = (x, u) -> x .+ Δt .* f_ode(x .+ (Δt / 2) .* f_ode(x, u), u)
dsys = MS.ConstrainedBlackBoxControlDiscreteSystem(f, 4, 1, _X_, _U_)

# Both angles periodic.
wrap = UT.get_periodic_wrapper(SVector(1, 2), SVector(2π, 2π); start = SVector(-π, -π))
T_split = UT.set_in_period(_T_, SVector(1, 2), SVector(2π, 2π), SVector(-π, -π))
discrete_problem =
    PR.OptimalControlProblem(dsys, _I_, T_split, nothing, nothing, PR.Infinity(), nothing)

x0 = SVector(0.0, 0.0, 0.0, 0.0)
nstep = 60      # his solution flips in 2.55 s — short & aggressive suits CEM
u_clamp = 7.5   # reserve 1.0 vs the plant's ±8.5

function rollout_inputs(us)
    xs = [x0]
    for u in us
        push!(xs, f(xs[end], u))
    end
    return ST.Trajectory(xs; inputs = collect(us))
end

rT = 0.9 .* [25π / 180, 25π / 180, 5.0, 5.0]
struct DPTerminalPull <: AB.AbstractCostTerm
    w::Float64
end
function AB.cost_final(t::DPTerminalPull, acc, xT)
    xw = wrap(xT)
    dθ1 = rem(xw[1] - π, 2π, RoundNearest)
    dθ2 = rem(xw[2] - π, 2π, RoundNearest)
    return acc +
           t.w * ((dθ1 / rT[1])^2 + (dθ2 / rT[2])^2 + (xw[3] / rT[3])^2 + (xw[4] / rT[4])^2)
end

# Up-up: total height of BOTH links (2 + cos θ₁ + cos θ₂, zero at the handstand).
struct HeightShapingCost <: AB.AbstractCostTerm
    w::Float64
end
AB.cost_step(t::HeightShapingCost, acc, x, u, k) = acc + t.w * (2 + cos(x[1]) + cos(x[2]))

# No pace cost: the handstand NEEDS the violent flip (ω domain ±6, target ±5).
score = AB.CompositeCost(
    AB.ReachObjectiveCost(T_split; wrap = wrap),
    DPTerminalPull(300.0),
    HeightShapingCost(5.0),
    AB.InputEffortCost(0.001),
    AB.InputSmoothnessCost(; w_du = 0.02, w_ddu = 0.005),
    AB.DomainPenaltyCost(_X_; wrap = wrap),
)

import Dionysos.Optim.Abstraction: rollout_cost
function scan_seeds()
    best_seed, best_cost, best_cfg = nothing, Inf, nothing
    for A in (2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
        ω in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0),
        φ in (0.0, π / 2, π)

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
best_seed, best_seed_cost, best_cfg = scan_seeds()
println("seed scan best: cost=$(round(best_seed_cost; digits=1))  (A, ω, φ)=$(best_cfg)")

mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = best_seed,
    nstep = nstep,
    nsamples = 3000,
    niter = 120,
    anneal = 0.99,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(1.5)),
    project_input = u -> SVector(clamp(u[1], -u_clamp, u_clamp)),
    cost = score,
    wrap_state = (p, x) -> wrap(x),
    update_rule = :cem,
    elite_frac = 0.05,
    antithetic = true,
    stop_on_success = false,
)
AB.set_problem!(mppi, discrete_problem)
gen_time = @elapsed AB.generate!(mppi)
traj = AB.get_trajectory(mppi)
xe = wrap(collect(ST.states(traj))[end])
println(
    "  $(round(gen_time; digits = 1)) s, success = $(AB.get_success(mppi)), ",
    "steps = $(length(ST.inputs(traj)))",
)
println("  endpoint (wrapped) = $(round.(collect(xe); digits = 3))")
AB.get_success(mppi) || error("generation failed — retune before certifying")

# Certification
lift = function (tr)
    lifted = ST.unwrap_trajectory(tr, (1, 2), (2π, 2π))
    xN = collect(ST.states(lifted))[end]
    s1 = 2π * round((xN[1] - π) / (2π))
    s2 = 2π * round((xN[2] - π) / (2π))
    (s1 == 0.0 && s2 == 0.0) && return lifted
    xs = [SVector(x[1] - s1, x[2] - s2, x[3], x[4]) for x in ST.states(lifted)]
    return ST.Trajectory(xs; inputs = collect(ST.inputs(lifted)))
end

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
    ST.format_input_set(_U_),
    ST.format_noise_set(Wset),
)

t = [0.15, 0.15, 0.75, 0.75]
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{4}(LazySets.low(H) ./ t),
    high = SVector{4}(LazySets.high(H) ./ t),
)
lifted_ref = lift(traj)
θ1s = [x[1] for x in ST.states(lifted_ref)]
θ2s = [x[2] for x in ST.states(lifted_ref)]
X_gate = LazySets.Hyperrectangle(;
    low = SVector(minimum(θ1s) - 1.0, minimum(θ2s) - 1.0, -6.0, -6.0),
    high = SVector(maximum(θ1s) + 1.0, maximum(θ2s) + 1.0, 6.0, 6.0),
)
println(
    "lifted corridors: θ₁ [$(round(minimum(θ1s); digits=2)), $(round(maximum(θ1s); digits=2))], ",
    "θ₂ [$(round(minimum(θ2s); digits=2)), $(round(maximum(θ2s); digits=2))]",
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
ztraj_fn(tr) = ST.Trajectory(
    [SVector{4}(collect(x) ./ t) for x in ST.states(tr)];
    inputs = collect(ST.inputs(tr)),
)

sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
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
    state_scaling = nothing,
    linearization_δx = [0.05, 0.05, 0.05, 0.05] ./ t,
    linearization_δu = [0.2],
    adaptive_boxes = adaptive_opts,
    objective = :maximin,
    remainder_model = :ball,
    domain_cap = true,
    check_state_domain = true,
)

println("— certify —")
flush(stdout)
bw = EB.BackwardCertifier(zprovider, sdp, back_opts)
EB.set_problem!(bw, zproblem)
EB.set_trajectory!(bw, ztraj_fn(lifted_ref))
cert_time = @elapsed EB.certify!(bw)

bres = EB.get_result(bw)
bres === nothing && error("no certification attempted")
println("  $(round(cert_time; digits = 1)) s, certified = $(bres.success)")
println(
    "  steps = $(length(ST.inputs(traj))), failed_k = $(bres.failed_k), ",
    "terminal ⊆ target: $(bres.terminal_contained), initial coverage: $(bres.initial_coverage)",
)
if !bres.success && !isempty(bres.steps)
    zx = collect(ST.states(ztraj_fn(lifted_ref)))
    for rec in bres.steps[1:min(3, length(bres.steps))]
        s = rec.summary
        xk = rec.k <= length(zx) ? round.(zx[rec.k] .* t; digits = 2) : "?"
        println(
            "  k=$(rec.k) [$(rec.status)] x=$(xk) ",
            "box_status=$(get(s, :adaptive_box_status, "?")) ",
            "cand=$(get(s, :candidate_statuses, "?"))",
        )
    end
end
# ------------------------------------------------------------
# Capture-phase certificate: re-certify the certified suffix standalone
# ------------------------------------------------------------

ztr = ztraj_fn(lifted_ref)
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

ells = cres.lmi_data.ellipsoids
vols = [π^2 / 2 * sqrt(LA.det(D * Matrix(LazySets.shape_matrix(E)) * D)) for E in ells]
phys_radii(E) = t .* sqrt.(LA.eigvals(LA.Symmetric(Matrix(LazySets.shape_matrix(E)))))
println("  capture depth profile (guaranteed recovery region vs time before target):")
for d in unique(clamp.([5, 10, 15, length(ells) - 1], 1, length(ells) - 1))
    E = ells[end - d]
    r = sort(phys_radii(E))
    println(
        "    $(lpad(round(d * Δt; digits = 2), 5)) s out: V = ",
        "$(round(vols[end - d]; sigdigits = 3)), semi-axes $(round.(r; sigdigits = 2))",
    )
end

zctrl = EB.get_controller(capture_bw)
println("— certified controller: FunnelController with $(length(zctrl.kappas)) steps —")

# Closed-loop validation from the meaningful-region depth.
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
        u = SVector{1}(clamp.(Matrix(κ.A) * collect(z) .+ collect(κ.c), -8.5, 8.5))
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
# Plots: the (θ₁, θ₂) corner run (the view of Florentin's animation) + the
# (θ₂, ω₂) swing plane, with funnel shadows
# ------------------------------------------------------------

denorm_shadow(E, i, j) = begin
    Q = D * Matrix(LazySets.shape_matrix(E)) * D
    c = t .* collect(LazySets.center(E))
    LazySets.Ellipsoid(c[[i, j]], Matrix(LA.Symmetric(Q[[i, j], [i, j]])))
end

xs_lift = collect(ST.states(lifted_ref))
cap_xs = xs_lift[capture_start:end]

fig1 = plot(;
    xlabel = "θ₁  [rad] (lifted)",
    ylabel = "θ₂  [rad] (lifted)",
    title = "Up-up handstand: certified capture corridor",
    legend = :topleft,
    size = (900, 650),
)
plot!(
    fig1,
    LazySets.Hyperrectangle(;
        low = SVector(π - 25π / 180, π - 25π / 180),
        high = SVector(π + 25π / 180, π + 25π / 180),
    );
    color = :green,
    alpha = 0.3,
    label = "target (θ₁, θ₂)",
)
for (i, E) in enumerate(ells)
    plot!(
        fig1,
        denorm_shadow(E, 1, 2);
        color = :steelblue,
        alpha = 0.25,
        linewidth = 1.0,
        linecolor = :steelblue,
        label = i == 1 ? "funnel (θ₁θ₂ shadow)" : "",
    )
end
plot!(
    fig1,
    [x[1] for x in xs_lift],
    [x[2] for x in xs_lift];
    color = :gray,
    linewidth = 1.5,
    label = "nominal (full flip)",
)
plot!(
    fig1,
    [x[1] for x in cap_xs],
    [x[2] for x in cap_xs];
    color = :black,
    linewidth = 2.5,
    label = "certified capture phase",
)

fig2 = plot(;
    xlabel = "θ₂  [rad] (lifted)",
    ylabel = "ω₂  [rad/s]",
    title = "swing plane",
    legend = :topleft,
)
plot!(
    fig2,
    LazySets.Hyperrectangle(;
        low = SVector(π - 25π / 180, -5.0),
        high = SVector(π + 25π / 180, 5.0),
    );
    color = :green,
    alpha = 0.3,
    label = "target (θ₂, ω₂)",
)
for E in ells
    plot!(
        fig2,
        denorm_shadow(E, 2, 4);
        color = :steelblue,
        alpha = 0.25,
        linewidth = 1.0,
        linecolor = :steelblue,
        label = "",
    )
end
plot!(
    fig2,
    [x[2] for x in xs_lift],
    [x[4] for x in xs_lift];
    color = :gray,
    linewidth = 1.5,
    label = "nominal",
)
plot!(
    fig2,
    [x[2] for x in cap_xs],
    [x[4] for x in cap_xs];
    color = :black,
    linewidth = 2.5,
    label = "capture phase",
)

final_fig = plot(fig1, fig2; layout = (1, 2), size = (1500, 650))
plot_path = joinpath(@__DIR__, "demo_double_pendulum_upup.png")
savefig(final_fig, plot_path)
println("— plot saved: $plot_path —")
