# # DC-DC converter: keeping a boost converter in a safe band
#
# | | |
# |:--|:--|
# | **System**        | 2-D [continuous](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem), [switched](https://blegat.github.io/HybridSystems.jl/stable/lib/types/#HybridSystems.AbstractHybridSystem) by the input |
# | **Specification** | [safety](@ref Dionysos.Problem.SafetyProblem) |
# | **Solver**        | [uniform grid abstraction](@ref Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer), twice |
#
# A boost DC-DC converter, widely studied as a hybrid control benchmark
# [girard2010approximately](@cite). The state is the inductor current and the capacitor voltage,
# $x = (i_L, v_C)$, and the only control is the position of a switch, which selects between two
# affine dynamics:
#
# ```math
# \dot{x} = A_p x + b, \qquad p \in \{1, 2\}.
# ```
#
# The goal is not to reach anything — it is to **never leave** a band around the operating
# point, using nothing but the switching signal.
#
# This page does one model twice. The plant and the specification are written once, and then
# handed to two different abstractions: a uniform grid, and a grid of ellipsoidal cells that
# exploits the converter's incremental stability. Nothing about the model changes between them,
# which is the point of putting every solver behind the same interface.

using StaticArrays, JuMP, Plots

# The first model supplies its dynamics as a Julia function and needs no symbolic backend; the
# hybrid model at the end writes them as `∂` expressions, which does.
using Symbolics, MathOptSymbolicAD

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction;

using Test     #src

# ## The model
#
# The circuit constants, and the two affine vector fields the switch selects between:
# ``\dot{x} = A_p x + b``, with ``p = 1`` for the closed switch and ``p = 2`` for the open one.

vs, xL, xC, r0, rL, rC = 1.0, 3.0, 70.0, 1.0, 0.05, 0.005

A1 = SMatrix{2, 2}(-rL / xL, 0.0, 0.0, -1.0 / (xC * (r0 + rC)))
A2 = SMatrix{2, 2}(
    -(rL + r0 * rC / (r0 + rC)) / xL,
    5.0 * r0 / ((r0 + rC) * xC),
    -r0 / ((r0 + rC) * xL * 5.0),
    -1.0 / (xC * (r0 + rC)),
)
b = SVector(vs / xL, 0.0)

dynamic = (x, u) -> u[1] == 1 ? A1 * x + b : A2 * x + b;

# The growth bound needs a *nonnegative* off-diagonal — it bounds how fast neighbouring
# trajectories separate, which is a magnitude — so mode 2 contributes `A2` with that entry
# made positive. Mode 1 is already diagonal and non-positive, so it is used as is.

A2_abs = SMatrix{2, 2}(
    -(rL + r0 * rC / (r0 + rC)) / xL,
    5.0 * r0 / ((r0 + rC) * xC),
    r0 / ((r0 + rC) * xL * 5.0),
    -1.0 / (xC * (r0 + rC)),
)

jacobian_bound = u -> u[1] == 1 ? A1 : A2_abs;

#-

x_low, x_upp = [1.15, 5.45], [1.55, 5.85]
safe = UT.box(SVector(x_low...), SVector(x_upp...))
initial = UT.box(SVector(1.19, 5.59), SVector(1.21, 5.61));

# `direct_model` rather than `Model`: the roles below are recorded against the optimizer, and a
# cached model has not handed it the variables yet.

model = direct_model(Dionysos.Optimizer());
@variable(model, x_low[i] <= x[i = 1:2] <= x_upp[i])
@variable(model, 1 <= u <= 2);

# The dynamics are a **Julia function**, not a JuMP expression: `Aₚ x + b` is selected by an
# `if` on the input, which no algebraic expression can express. Supplying `f` leaves the
# front-end nothing to infer from, so two things have to be said explicitly — which variables
# are states, and whether `f` is a vector field or a one-step map.

set_role!(x, Dionysos.STATE)
set_attribute(model, "dynamics", dynamic)
set_attribute(model, "time_domain", Dionysos.CONTINUOUS);

# Safety: start in a small box around the operating point, and never leave the band. `Always`
# rather than `∉` — staying inside is the *requirement*, so those states must remain
# representable for the synthesis to reason about leaving them.

@constraint(model, x in Start(initial))
@constraint(model, x in Always(safe));

# ## The abstraction
#
# The input grid is what makes the switch a two-valued signal: one cell wide, centred on 1,
# over the input bounds `[1, 2]`.

hx = SVector(2.0 / 4.0e3, 2.0 / 4.0e3)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), hx))
set_attribute(model, "input_grid", MP.GridFree(SVector(1), SVector(1)))
set_attribute(model, "jacobian_bound", jacobian_bound)
set_attribute(model, "time_step", 0.5)
set_attribute(model, "approx_mode", AB.UniformGridAbstraction.GROWTH)
set_attribute(model, "print_level", 0);

# ## Solving

optimize!(model);

# ## Results

termination_status(model)

#-

concrete_problem = get_attribute(model, "concrete_problem");
concrete_controller = get_attribute(model, "concrete_controller");
invariant_set = get_attribute(model, "invariant_set");

@test is_solved_and_feasible(model)     #src

# The certificate of a safety problem is the **maximal controlled-invariant set**: every state
# in it can be kept inside the band forever, and no state outside it can.

# ## Closed loop

trajectory = Dionysos.simulate(model, SVector(1.2, 5.6); nsteps = 150);

@test all(x ∈ safe for x in ST.states(trajectory))     #src

# ## Visualisation

fig = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(trajectory; ms = 2.0, color = :blue)

# The converter circuit alongside the state and the switching signal. Only the drawing is
# borrowed from the problem library — everything above is defined on this page.

include(
    joinpath(dirname(dirname(pathof(Dionysos))), "problems", "DCDC", "dcdc_converter.jl"),
);

# `frame_step` and `fps` trade smoothness against file size. Two steps per frame keeps the
# state moving a small distance between frames — at six it jumped far enough to read as
# stuttering — and `fps` is chosen so the run lasts about ten seconds.

anim = Dionysos.animate_trajectory_dashboard(
    DCDC.system_plot!(),
    trajectory;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.5,
    frame_step = 2,
    xlabel_state = "iL [A]",
    ylabel_state = "vC [V]",
    ylabel_input = "switch position",
);
gif(anim; fps = 8)

# ## The same problem, a different abstraction
#
# `concrete_problem` is a solver-independent object: a system plus a
# [`SafetyProblem`](@ref Dionysos.Problem.SafetyProblem). Handing it to another optimizer needs
# no change to the model at all — only a different discretisation.
#
# This converter is *incrementally stable*, which means trajectories from nearby states
# converge. That justifies covering the state space with ellipsoidal cells shaped by a
# Lyapunov matrix `P`, instead of squares: fewer, better-adapted cells for the same guarantee.

origin = SVector(0.0, 0.0)
η = (2 / 4.0) * 10^(-3)
ϵ = 0.1 * 0.01
P = SMatrix{2, 2}(1.0224, 0.0084, 0.0084, 1.0031)
state_grid = MP.GridEllipsoidalRectangular(origin, SVector(η, η), P / ϵ)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("input_grid"),
    MP.GridFree(SVector(1), SVector(1)),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.5)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
MOI.optimize!(optimizer);

#-

MOI.get(optimizer, MOI.TerminationStatus())

#-

abstract_system2 = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"));
controller2 = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"));

traj2 = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    controller2,
    SVector(1.2, 5.6),
    300,
);

@test all(x ∈ safe for x in ST.states(traj2))     #src

fig2 = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(traj2; ms = 2.0, color = :blue)

# ## The same plant, a different model
#
# Everything above treats the converter as **one continuous system whose input selects the
# dynamics**. But the plant is just as honestly described as a **hybrid automaton**: two modes,
# one per switch position, each carrying its own affine dynamics, with transitions between them.
#
# Neither reading is a workaround for the other — they are two ways of saying the same thing,
# and the front-end expresses both. What changes is not the plant but the vocabulary: `@mode`
# instead of an input dimension, and `add_transition!` instead of `if u == 1`.
#
# Note what is *absent*: no input variable. Here the switch is the only control, so no mode has a
# continuous input at all. Its input space is zero-dimensional — a space with exactly one point,
# which is the action "leave the switch where it is" and is what lets the state evolve inside a
# mode. Nothing needs to be declared for it.

hybrid = Model(Dionysos.Optimizer);
@variable(hybrid, 1.15 <= xh[i = 1:2] <= 1.55)
set_lower_bound(xh[2], 5.45)
set_upper_bound(xh[2], 5.85)

@mode(hybrid, closed)
@mode(hybrid, opened);

# One mode per switch position, each with its own affine dynamics.

@constraint(closed, ∂(xh[1]) == A1[1, 1] * xh[1] + A1[1, 2] * xh[2] + b[1])
@constraint(closed, ∂(xh[2]) == A1[2, 1] * xh[1] + A1[2, 2] * xh[2] + b[2])
@constraint(opened, ∂(xh[1]) == A2[1, 1] * xh[1] + A2[1, 2] * xh[2] + b[1])
@constraint(opened, ∂(xh[2]) == A2[2, 1] * xh[1] + A2[2, 2] * xh[2] + b[2]);

# A guard spanning the whole safe set means the switch is always available — which is exactly
# what a converter's switch is: the controller decides, the state does not gate it.

add_transition!(hybrid, closed => opened) do t
    return @constraint(t, xh in Guard(safe))
end
add_transition!(hybrid, opened => closed) do t
    return @constraint(t, xh in Guard(safe))
end

@constraint(closed, xh in Always(safe))
@constraint(opened, xh in Always(safe));

# The discretisation is per mode, as always for a hybrid model. This one is deliberately four
# times coarser than the grid used above: the point of this section is the *encoding*, and the
# fine grid is already exercised by the first model.

for (md, Ab) in ((closed, A1), (opened, A2_abs))
    set_attribute(md, "state_grid", MP.GridFree(SVector(0.0, 0.0), 4 .* hx))
    set_attribute(md, "time_step", 0.5)
    set_attribute(md, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(md, "jacobian_bound", u -> Ab)
    set_attribute(md, "print_level", 0)
end

# `print_level` on a mode silences that mode's sub-solver; the hybrid solver above them has
# its own, set on the model.

set_attribute(hybrid, "print_level", 0)

optimize!(hybrid);

#-

termination_status(hybrid)

#-

@test is_solved_and_feasible(hybrid)     #src

# The augmented state of a hybrid model is `(x, mode)`, so the closed loop starts from a state
# *and* a switch position, and carries the mode the controller chose at each step alongside the
# usual state and input channels.

hybrid_traj = Dionysos.simulate(hybrid, (SVector(1.2, 5.6), 1); nsteps = 200);

@test all(x ∈ safe for x in ST.states(hybrid_traj))     #src

# Projected back onto the state plane, the hybrid closed loop is directly comparable with the
# two above: same band, same operating point, a third route to staying inside it.

fig3 = plot(; aspect_ratio = :equal)
plot!(concrete_problem)
plot!(hybrid_traj; ms = 2.0, color = :blue)

# And the same dashboard as the first model. One translation is needed first: the *inputs* of a
# hybrid closed loop are symbolic — an empty vector for "hold the switch", a `"SWITCH 1 -> 2"`
# marker when it flips — so there is no numeric channel for the input panel to draw. What the
# panel should show is the switch *position*, and that is the `modes` channel. Re-expressing the
# mode as a one-dimensional input gives back exactly the signal the continuous model plotted.

hybrid_switch = ST.Trajectory(
    ST.states(hybrid_traj);
    inputs = [SVector(Float64(m)) for m in ST.modes(hybrid_traj)[1:(end - 1)]],
);

anim_hybrid = Dionysos.animate_trajectory_dashboard(
    DCDC.system_plot!(),
    hybrid_switch;
    xdims = (1, 2),
    udims = (1,),
    Δt = 0.5,
    frame_step = 2,
    xlabel_state = "iL [A]",
    ylabel_state = "vC [V]",
    ylabel_input = "switch position",
);
gif(anim_hybrid; fps = 8)

# ## References
# 1. A. Girard, G. Pola and P. Tabuada, "Approximately Bisimilar Symbolic Models for Incrementally Stable Switched Systems," in IEEE Transactions on Automatic Control, vol. 55, no. 1, pp. 116-126, Jan. 2010.
# 2. S. Mouelhi, A. Girard, and G. Gössler. “CoSyMA: a tool for controller synthesis using multi-scale abstractions”. In: HSCC. ACM. 2013, pp. 83–88.
# 3. A. Girard. “Controller synthesis for safety and reachability via approximate bisimulation”. In: Automatica 48.5 (2012), pp. 947–953.
