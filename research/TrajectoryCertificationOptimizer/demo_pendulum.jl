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
# Status (2026-08-11, Clarabel, fixed rng): FULLY CERTIFIED. Seed 34 states in
# ~23 s; the driver loop certifies in round 1 (~56 s): CEM-MPPI finds the 48-step
# swing-up and the backward chain certifies ALL 48 steps with the terminal gate
# passing, on the plant's true ±4.5 input set — the P0 baseline failed at k≈19
# against an unsound ±10.5 input set. The demo ends with a 48-step
# ST.FunnelController. Honest open item, reported by the gates: initial coverage
# ≫ 1 (the entry funnel does not cover the full ±10°×±0.5 initial set — gap D;
# the forward direction / entry enlargement is the designed answer), and the
# forward↔backward handoff does not yet nest on this config.

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
u_plan = 0.75 * u_max            # input reserve: leave headroom for the feedback

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
    nstep = 60,
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
# The driver runs generation; between generation and certification it applies the
# periodic lift, shifted so the endpoint lands in the target's θ-range.
prepare = function (traj)
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

# The unwrapped cover has no θ seam and the objective has no obstacle; θ spans the
# full period, so the domain gate would only re-check |ω| ≤ 7 while flagging the
# lifted θ — disable it here and rely on the generous ω margin.
adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    true,
    [0.05, 0.10],
    [0.005, 0.005],
    [1.8, 2.5],
    [0.25],
    [0.01],
    [4.5],
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
    maxδx = 1.5,
    maxδu = 3.0,
    λ = 0.001,
    terminal_shape = nothing,            # default: inscribed ellipsoid of the target
    terminal_shrink = 0.85,
    state_scaling = LA.diagm([0.85 * 15.0 * π / 180.0, 0.25]),
    linearization_δx = [0.2, 0.4],
    linearization_δu = [1.0],
    adaptive_boxes = adaptive_opts,
    objective = :maximin,
    check_state_domain = false,
)

println("— generate ⇄ certify loop (retry ladder, up to 5 rounds) —")
bw = EB.BackwardCertifier(provider, sdp, back_opts)
driver = AB.TrajectoryCertificationOptimizer.Optimizer(mppi, bw)
MOI.set(driver, MOI.RawOptimizerAttribute("concrete_problem"), discrete_problem)
MOI.set(driver, MOI.RawOptimizerAttribute("max_rounds"), 5)
MOI.set(driver, MOI.RawOptimizerAttribute("prepare_trajectory"), prepare)
loop_time = @elapsed MOI.optimize!(driver)

rounds = MOI.get(driver, MOI.RawOptimizerAttribute("rounds"))
loop_success = MOI.get(driver, MOI.RawOptimizerAttribute("success"))
traj = MOI.get(driver, MOI.RawOptimizerAttribute("trajectory"))
bres = EB.get_result(bw)
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
lifted = prepare(traj)

println("— forward certification + handoff —")
fwd_opts = EB.ForwardOptions(;
    target_mode = :free,
    q_min = 1e-8,
    q_max = 1e2,
    maxδu = 3.0,
    linearization_δu = [1.0],
    λ = 0.001,
    # A realistic entry tube: the full circumscribed initial set (radius ≈ 0.7 in
    # ω) is not one-step certifiable with this input authority — the handoff demo
    # starts from a tighter but non-trivial entry neighbourhood.
    entry_shape = LA.diagm([0.05^2, 0.1^2]),
    check_state_domain = false,
)
fw = EB.ForwardCertifier(provider, sdp, fwd_opts)
bw2 = EB.BackwardCertifier(provider, sdp, back_opts)
ho_time = @elapsed handoff = EB.bidirectional_certify!(fw, bw2, problem, lifted)
fres = handoff.forward_result
println(
    "  $(round(ho_time; digits = 1)) s, forward steps certified: ",
    fres.failed_k === nothing ? length(ST.inputs(lifted)) : fres.failed_k - 1,
    ", handoff = $(handoff.success) at k = $(handoff.k_handoff)",
)

if handoff.success || loop_success
    ctrl = loop_success ? AB.get_controller(bw) : handoff.controller
    println(
        "— certified controller: $(typeof(ctrl).name.name) with ",
        "$(length(ctrl.kappas)) steps —",
    )
end
