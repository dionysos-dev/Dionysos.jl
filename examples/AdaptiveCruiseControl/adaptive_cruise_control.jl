# Adaptive cruise control — one abstraction, two specifications.
#
# The plant is a longitudinal point mass following a car travelling at a constant speed; the
# state is `(gap, ego speed)` and the input is an acceleration. The abstraction is built once
# and reused: switching the control task does not rebuild it.

using Dionysos

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "AdaptiveCruiseControl",
        "adaptive_cruise_control.jl",
    ),
)

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const ACC = AdaptiveCruiseControl

using Plots
import LazySets
import MathOptInterface as MOI
import StaticArrays: SVector

params = ACC.Params()

# Well behind the lead and travelling at about its speed, so the gap neither opens nor closes
# on its own: whatever the ego does to it, the controller chose to do.
x0 = SVector(90.0, 13.0)
_I_ = LazySets.Hyperrectangle(; low = SVector(89.0, 12.8), high = SVector(91.0, 13.2))

# `time_step ⋅ input step = 0.2 = speed step`, so every input translates the speed axis onto
# itself by a whole number of cells. The gap axis is not aligned, and the resistance term would
# break exactness anyway. Halving both state steps costs roughly 4× everywhere.
Δt = 0.5
state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 0.2))
input_grid = MP.GridFree(SVector(0.0), SVector(0.4))

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), ACC.jacobian_bound(params))
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

# ------------------------------------------------------------
# 1. Safety: never tailgate
# ------------------------------------------------------------

safety_problem = ACC.safety_problem(; params = params)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), safety_problem)
MOI.optimize!(optimizer)
println("safety: ", MOI.get(optimizer, MOI.TerminationStatus()))

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
XMapping = SY.get_state_mapping(abstract_system)
invariant_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("invariant_set"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

# The computed invariant set against the closed form it can be checked with. Braking without
# the resistance force needs more gap than is really necessary, so that curve is an inner
# bound; braking with it gives an outer one. An invariant set reaching below the outer curve
# would mean the abstraction is unsound.
speeds = range(params.v_min, params.v_max; length = 200)
a_outer = ACC.max_deceleration(params, params.v_max)

fig = plot(; xlabel = "gap [m]", ylabel = "ego speed [m/s]", legend = :bottomright)
plot!(fig, safety_problem)
plot!(
    fig,
    (invariant_set, XMapping);
    color = :blue,
    linecolor = :blue,
    label = "Invariant set",
)
plot!(
    fig,
    [ACC.min_safe_gap(params, v) for v in speeds],
    speeds;
    lw = 2,
    color = :black,
    label = "inner bound",
)
plot!(
    fig,
    [ACC.min_safe_gap(params, v; a_brake = a_outer) for v in speeds],
    speeds;
    lw = 2,
    ls = :dash,
    color = :black,
    label = "outer bound",
)
display(fig)

# ------------------------------------------------------------
# 2. Reach and stay: cruise at the set-point, or follow the lead
# ------------------------------------------------------------

# The disjunction is the specification, not a convenience. With the set-point at 24 m/s behind a
# car doing 13.89, the cruise half is unreachable — every gap in the envelope closes in finite
# time — so the run settles in the follow half. Handing the cruise band over on its own returns
# `LOCALLY_INFEASIBLE`, which is the benchmark proving something true rather than failing.
#
# `gap_high` is what makes the follow half a *following distance* rather than a floor. Left
# open, sitting 90 m back at the lead's own speed already satisfies "matched speed at a safe
# gap", and the ego never closes it. Capped, the ego has to overspeed, catch up, and then
# settle — the behaviour worth watching.
settle = UT.set_union([ACC.cruise_set(params), ACC.follow_set(params; gap_high = 45.0)])
ras_problem = ACC.reach_and_stay_problem(; params = params, _I_ = _I_, _T_ = settle)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), ras_problem)
MOI.optimize!(optimizer)
println(
    "reach-and-stay, set-point or follow: ",
    MOI.get(optimizer, MOI.TerminationStatus()),
)

ras_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
ras_traj = ST.get_closed_loop_trajectory(discrete_time_system, ras_controller, x0, 160)
let xs = collect(ST.states(ras_traj))
    println(
        "  gap ",
        round(xs[1][1]; digits = 1),
        " m -> ",
        round(xs[end][1]; digits = 1),
        " m (min ",
        round(minimum(x[1] for x in xs); digits = 1),
        "), top speed ",
        round(maximum(x[2] for x in xs); digits = 1),
        " m/s",
    )
end

winning_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("winning_set"))

fig = plot(; xlabel = "gap [m]", ylabel = "ego speed [m/s]", legend = :bottomright)
plot!(fig, ras_problem)
plot!(
    fig,
    (winning_set, XMapping);
    color = :yellow,
    linecolor = :yellow,
    label = "Winning set",
)
plot!(fig, ras_traj; ms = 2.0, color = :blue)
display(fig)

# There is no reachability-only run here on purpose. Reach-and-stay already has to reach, and
# it keeps going once it arrives; dropping the stay half would only stop the trajectory at
# first entry and show less.

# ------------------------------------------------------------
# The road view
# ------------------------------------------------------------

anim = Dionysos.animate_trajectory_dashboard(
    ACC.system_plot!(; params = params),
    ras_traj;
    xdims = (1, 2),
    udims = (1,),
    Δt = Δt,
    xlabel_state = "gap [m]",
    ylabel_state = "ego speed [m/s]",
    xlabel_input = "time [s]",
    ylabel_input = "acceleration [m/s²]",
)
gif(anim; fps = 4)
