# Certified articulated-vehicle maneuver — the P6 vehicle demo (plan.md §8.2).
#
# 4-D tractor-trailer (x, y, θ, ϕ) driving from the lower-left corner to the
# upper-right target. The pipeline runs without an abstraction seed: a constant-
# input rollout seeds CEM-MPPI, and certification happens backward-only in
# globally normalized coordinates with the state-domain gate ACTIVE (every funnel
# ellipsoid must stay inside X). Certified against the plant's true input set.
#
#     julia --project=bench research/TrajectoryCertificationOptimizer/demo_vehicle.jl
#
# Status (2026-08-11, Clarabel, fixed rng): generation SOLVED (110 steps reaching
# the target with |θ|, |ϕ| ≤ 0.2 — needed the terminal-ellipsoid shaping and a
# horizon long enough for the hitch angle's ~14-step relaxation), certification
# PARTIAL: the backward chain certifies steps 110→41 (69 steps, terminal gate
# passing, domain gate active) and fails mid-turn at k = 41 with
# `lmi_infeasible_at_max_box` — peak curvature + steering near the planning bound
# leaves too little feedback headroom. Next levers (plan.md): more steering
# reserve (u_plan δ 0.5 → 0.4), smaller Δt through the turn, the ball remainder
# model (16 vertex blocks per 4-D step), and prefix re-planning into E_42.

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

params = AV.Params()
problem = AV.problem(; params = params)
Δt = 0.5

base = PR.discretize_problem(problem, Δt; num_substeps = 5)
f = MS.mapping(base.system)

# ------------------------------------------------------------
# 1) Seed: constant-input rollout (no abstraction needed for this task)
# ------------------------------------------------------------

x0 = SVector(-9.5, -9.5, 0.0, 0.0)
# ≥ 26.9 m diagonal at v ≤ 0.85·Δt per step ⟹ ≥ 64 straight-line steps; the curve
# and the straight final approach that relaxes the hitch angle need margin (the
# trailer's φ has a ~14-step relaxation constant at this speed).
nstep = 110
seed_traj = begin
    u0 = SVector(0.85, 0.05)
    xs = [x0]
    for _ in 1:nstep
        push!(xs, f(xs[end], u0))
    end
    ST.Trajectory(xs; inputs = fill(u0, nstep))
end

# ------------------------------------------------------------
# 2) CEM-MPPI
# ------------------------------------------------------------

println("— MPPI —")
u_plan = SVector(0.85, 0.5)      # input reserve vs U = [-1,1]×[-0.6,0.6]

# Terminal shaping through the target's inscribed ellipsoid: the reach distance is
# position-dominated, but the target also demands θ, ϕ ∈ [−0.2, 0.2] — arriving
# along the diagonal means θ ≈ π/4 unless the endpoint alignment is *rewarded*.
E_T = LazySets.Ellipsoid(
    collect(LazySets.center(problem.target_set)),
    Matrix(LA.Diagonal([0.5, 0.5, 0.12, 0.12] .^ 2)),
)

cost = AB.CompositeCost(
    AB.ReachObjectiveCost(problem.target_set),
    AB.TerminalEllipsoidCost(E_T; w_outside = 1e5, w_center = 1e4),
    AB.InputEffortCost(0.01),
    AB.InputSmoothnessCost(; w_du = 0.5, w_ddu = 0.1),
    AB.DomainPenaltyCost(problem.system.X),
)

mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.MersenneTwister(1),
    seed_trajectory = seed_traj,
    nstep = nstep,
    nsamples = 1000,
    niter = 70,
    noise = AB.MPPITrajectoryGenerator.GaussianMPPINoise(SVector(0.2, 0.15)),
    project_input = u ->
        SVector(clamp(u[1], -u_plan[1], u_plan[1]), clamp(u[2], -u_plan[2], u_plan[2])),
    cost = cost,
    update_rule = :cem,
    elite_frac = 0.05,
    antithetic = true,
)

# ------------------------------------------------------------
# 3) Normalized symbolic provider (plan.md §4.3)
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
    ST.format_input_set(problem.system.U),   # the plant's true input set
    ST.format_noise_set(Wset),
)

t = [2.0, 2.0, 0.35, 0.35]                   # characteristic scales (m, m, rad, rad)
zprovider = ST.normalized_symbolic_provider(provider, t)
D = LA.Diagonal(t)

zbox(H) = LazySets.Hyperrectangle(;
    low = SVector{4}(LazySets.low(H) ./ t),
    high = SVector{4}(LazySets.high(H) ./ t),
)
# The z-frame problem carries the scaled domain so the state-domain gate checks
# every funnel ellipsoid against X in the same coordinates it was certified in.
zsys = AV.system(zbox(problem.system.X); params = params)
zproblem = PR.OptimalControlProblem(
    zsys,
    zbox(problem.initial_set),
    zbox(problem.target_set),
    nothing,
    nothing,
    base.time,
    nothing,
)

ztraj(traj) = ST.Trajectory(
    [SVector{4}(collect(x) ./ t) for x in ST.states(traj)];
    inputs = collect(ST.inputs(traj)),
)

# ------------------------------------------------------------
# 4) Backward certification with the domain gate ON
# ------------------------------------------------------------

sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    true,
    [0.05, 0.05, 0.05, 0.05] ./ t,
    [0.005, 0.005, 0.005, 0.005] ./ t,
    [3.0, 3.0, 1.0, 1.0] ./ t,
    [0.2, 0.15],
    [0.02, 0.015],
    [1.0, 0.6],
    1.5,
    1.05,
    30,
    1e-8,
    false,
    [1.0],
    :first_consistent,
    true,
)
back_opts = EB.ChainOptions(;
    maxδx = 8.0,
    maxδu = 1.0,
    λ = 0.001,
    terminal_shape = nothing,
    terminal_shrink = 0.9,
    state_scaling = nothing,                 # exact: dynamics already normalized
    linearization_δx = [0.05, 0.05, 0.05, 0.05] ./ t,
    linearization_δu = [0.3, 0.2],
    adaptive_boxes = adaptive_opts,
    objective = :maximin,
    check_state_domain = true,
)

println("— generate ⇄ certify loop (retry ladder, up to 3 rounds) —")
bw = EB.BackwardCertifier(zprovider, sdp, back_opts)
driver = AB.TrajectoryCertificationOptimizer.Optimizer(mppi, bw)
MOI.set(driver, MOI.RawOptimizerAttribute("concrete_problem"), base)
MOI.set(driver, MOI.RawOptimizerAttribute("certifier_problem"), zproblem)
MOI.set(driver, MOI.RawOptimizerAttribute("max_rounds"), 3)
MOI.set(driver, MOI.RawOptimizerAttribute("prepare_trajectory"), ztraj)
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
    "$(bres.terminal_contained), initial coverage: $(bres.initial_coverage), ",
    "domain gate ran: $(bres.state_domain_checked)",
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

if loop_success
    zctrl = AB.get_controller(bw)
    println("— certified controller: FunnelController with $(length(zctrl.kappas)) steps —")
end

# ------------------------------------------------------------
# 5) Plot in the (x, y) plane: trajectory + funnel shadows
# ------------------------------------------------------------

fig = plot(;
    xlabel = "x  [m]",
    ylabel = "y  [m]",
    title = "Certified articulated-vehicle maneuver",
    legend = :topleft,
    size = (900, 700),
    aspect_ratio = :equal,
)

box2(H) = LazySets.Hyperrectangle(;
    low = SVector(LazySets.low(H, 1), LazySets.low(H, 2)),
    high = SVector(LazySets.high(H, 1), LazySets.high(H, 2)),
)
plot!(fig, box2(problem.initial_set); color = :gray, alpha = 0.5, label = "initial set")
plot!(fig, box2(problem.target_set); color = :green, alpha = 0.35, label = "target set")

# xy-shadow of a 4-D funnel ellipsoid: the principal 2×2 block of the
# de-normalized shape matrix.
function funnel_shadow(E)
    Q = LA.Diagonal(t) * Matrix(LazySets.shape_matrix(E)) * LA.Diagonal(t)
    c = t .* collect(LazySets.center(E))
    return LazySets.Ellipsoid(c[1:2], LA.Symmetric(Q[1:2, 1:2]) |> Matrix)
end

if bres !== nothing
    for (i, E) in enumerate(bres.lmi_data.ellipsoids)
        plot!(
            fig,
            funnel_shadow(E);
            color = :steelblue,
            alpha = 0.2,
            linewidth = 0,
            label = i == 1 ? "funnel (xy shadow)" : "",
        )
    end
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

plot_path = joinpath(@__DIR__, "demo_vehicle.png")
savefig(fig, plot_path)
println("— plot saved: $plot_path —")
