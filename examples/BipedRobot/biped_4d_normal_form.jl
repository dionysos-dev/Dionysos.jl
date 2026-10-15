# # Underactuated four-link biped in normal form: one footstep over an obstacle
#
# Same scenario as `biped_4d_underactuated.jl` (a footstep over a step at the stance foot, with
# the three motor joints commanded in velocity and the ankle `q₁` unactuated), but with a
# different model and a different abstraction:
#
# - **model**: the multibody system is built *symbolically* by `RP.symbolic_model()`
#   (`problems/BipedRobot/4D_model_underactuated/robot_problem.jl`, on top of `Multibody2D.jl`; no
#   RigidBodyDynamics in the dynamics), brought to normal form `x = (q, v_a, p_u)`, and reduced to the
#   velocity-controlled design model in which the unactuated ankle is described by its
#   **generalized momentum** `p₁ = ∂L/∂q̇₁` instead of its velocity (`RP.reduced_dynamics`):
#
#       x = [p₁, q₁, q₂, q₃, q₄],   u = [v₂, v₃, v₄],
#
#   i.e. `q̇₁ = (p₁ - D₁ₐ v_a) / D₁₁`, `ṗ₁ = ∂L/∂q₁`, `(q̇₂, q̇₃, q̇₄) = v_a`, with the commands
#   held constant over a sampling period. `p₁` takes the slot of `q̇₁` in the original state, so the
#   foot-height constraint, the obstacle carving, the target set and the plots of
#   `robot_problem.jl` apply unchanged;
# - **abstraction**: a certified over-approximation with the **local growth bound** of
#   `LocalGrowthAbstraction.jl` (a Jacobian bound per cell and input, with interval arithmetic and
#   a verified flow enclosure).
#   The abstraction is nondeterministic, so the controller is synthesized against the worst case.
#
# The constraints are those of the original example: the swing-foot height is bounded below by the
# ground (cells where the foot is surely below `z = 0` are removed) and the triangular obstacle at
# the origin is carved out, which forces the foot to lift.
#
# Conventions: `Multibody2D` and `robot_problem.jl` both measure angles counter-clockwise, so the
# symbolic model is used as is. The forward kinematics of the carving and of the target set are
# compiled from the same symbolic model (`RP.Geometry`), so the dynamics and the constraints cannot
# disagree. The mechanism is the one of `four_link_biped.urdf` (thin rods `m L² / 12`, point-mass
# motors), the same as in `biped_4d_underactuated.jl`.
#
# Run it with (an environment with Dionysos, Symbolics, IntervalArithmetic, Plots and the other
# packages of `problems/BipedRobot`; RigidBodyDynamics is not needed):
#
# ```
# julia -t auto --project=problems/BipedRobot examples/BipedRobot/biped_4d_normal_form.jl
# ```

using StaticArrays
using MathematicalSystems
using Dionysos
import MathOptInterface as MOI
import LazySets

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

const PROBLEM_DIR = joinpath(
    dirname(dirname(pathof(Dionysos))),
    "problems",
    "BipedRobot",
    "4D_model_underactuated",
)

# Mechanism, geometry, obstacle carving, target set and plots (state `[ankle, q₁, q₂, q₃, q₄]`).
include(joinpath(PROBLEM_DIR, "robot_problem.jl"))
import .UnderactuatedRobotProblem as RP

# Local growth bound.
include(joinpath(PROBLEM_DIR, "LocalGrowthAbstraction.jl"))
import .LocalGrowthAbstraction as LGA

# ------------------------------------------------------------
# 1) Symbolic model → numerical dynamics
# ------------------------------------------------------------

values = RP.default_parameters()   # lengths, masses and gravity of four_link_biped.urdf

@time model = RP.symbolic_model()
println("Symbolic normal form: ", model.nf)
reduced = RP.reduced_dynamics(model)       # x = (p₁, q₁…q₄), u = (v₂, v₃, v₄), parameters still symbolic
f_sym = RP.substitute_parameters(model, reduced.f, values)   # now only x and u are free
x_sym, u_sym = reduced.x, reduced.u

# The forward kinematics of the carving and of the target set, compiled from the same model.
geometry = RP.Geometry(model, values)

# ------------------------------------------------------------
# 2) Scenario: the `step` of the fully actuated example
# ------------------------------------------------------------
# An 8 cm × 3 cm step at the stance foot that the swing foot must clear, and a
# foothold 20 cm in front of the stance foot.

_triangle(cx, hw, h) =
    LazySets.VPolygon([SVector(cx - hw, 0.0), SVector(cx, h), SVector(cx + hw, 0.0)])

scenario_name = :step_normal_form
obstacle = @isdefined(OBSTACLE) ? OBSTACLE : [_triangle(0.0, 0.04, 0.03)]
foothold = SVector(0.2, 0.0)

# Set `DISC` before including this file to try another resolution. The motor joints move by
# exactly one cell (`dx`) per step at speed `du = dx / tstep`; `dq̇1` is here the step of the
# ankle *momentum* `p₁` (`D₁₁ ≈ 0.14 kg m²`, so 0.02 corresponds to ≈ 0.15 rad/s). The step is
# shorter than in the original example (0.05 s instead of 0.1 s): the local growth bound needs
# `tstep · ‖∂f/∂x‖ ≲ 1` (here `‖∂f/∂x‖ ≈ 9`) for the flow enclosure to exist.
disc =
    @isdefined(DISC) ? DISC :
    RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.02, tstep = 0.05)

# State box `[p₁, q₁, q₂, q₃, q₄]`: narrower than that of the original example, to keep the
# abstraction affordable at the fine ankle resolution above (287 793 cells before carving, 176 225
# after, i.e. 4.8 M cell/input pairs: about two minutes on 10 threads). `|p₁| ≤ 0.2 N m s` covers `|q̇₁| ≲ 1.2 rad/s`
# (`p₁ ≈ D₁₁ q̇₁ + D₁ₐ v_a`, `D₁₁ ≈ 0.14`). Angles are counter-clockwise: the swing hip swings
# forwards for `q₃ > 0` (up to 1.0 here, to leave room for the start posture below), the swing
# knee bends backwards for `q₄ < 0`, and the robot falls forwards, towards negative `q₁`
# (hence the room on that side). Set `X_BOUNDS = (x_lb, x_ub)` to change it, e.g. to the wider
# `(SVector(-0.3, -0.4, -0.8, -0.8, -1.2), SVector(0.3, 0.4, 0.8, 1.2, 0.2))` (3× more cells).
x_lb, x_ub =
    @isdefined(X_BOUNDS) ? X_BOUNDS :
    (SVector(-0.2, -0.3, -0.5, -0.8, -1.0), SVector(0.2, 0.2, 0.5, 1.0, 0.2))

# The footstep starts from the configuration `q0 = [q₁, q₂, q₃, q₄]` and ends on the foothold, with
# **no velocity imposed**: the controller must work from every value of the first state coordinate
# (the ankle momentum `p₁`, which replaces the ankle velocity) in `INITIAL_VELOCITY` and may end with
# any value in `TARGET_VELOCITY`. Both are `(vlo, vhi)` ranges on `p₁` and default to its whole
# interval of definition (`nothing`). The swing foot has no ground reaction force, so requiring rest
# is needlessly conservative. Set `Q0`, `INITIAL_VELOCITY` or `TARGET_VELOCITY` before including this
# file to change them, e.g. `INITIAL_VELOCITY = (-0.05, 0.05)`.
#
# The swing foot has to start **outside** the carved neighbourhood of the obstacle: the carving is
# sound, so it removes every cell whose swing foot is within the Lipschitz margin (printed below) of
# the obstacle. The default has the swing leg straight and 0.6 rad behind the hip, which puts the
# foot at x ≈ -0.17 m, 13 cm from the obstacle's edge.
q0 = @isdefined(Q0) ? Q0 : SVector(-0.1, 0.0, -0.5, 0.0)
initial_velocity = @isdefined(INITIAL_VELOCITY) ? INITIAL_VELOCITY : (0.5, 0.5)
target_velocity = @isdefined(TARGET_VELOCITY) ? TARGET_VELOCITY : nothing

println("Carving the obstacle and the ground out of the state domain…")
X_box = LazySets.Hyperrectangle(; low = x_lb, high = x_ub)
U_box = LazySets.Hyperrectangle(;
    low = SVector(-disc.u_max, -disc.u_max, -disc.u_max),
    high = SVector(disc.u_max, disc.u_max, disc.u_max),
)
@time removed = RP.infeasible_cells(X_box, disc.state_grid, obstacle, geometry)
margin = RP.swing_foot_deviation_bound(geometry, MP.get_h(disc.state_grid) ./ 2)
println("Lipschitz carving margin: ", round(100 * margin; digits = 1), " cm")
X = UT.set_minus(X_box, removed)

# ------------------------------------------------------------
# 3) Abstraction with local growth bounds
# ------------------------------------------------------------

@time lgb = LGA.LocalGrowthBound(f_sym, x_sym, u_sym, X, U_box, disc.tstep)
concrete_system = ST.get_system(lgb)   # discrete-time (RK4) system on the carved domain

# Concrete start of the closed-loop run below: the configuration `q0` with the first coordinate in the
# middle of the initial range, snapped to the grid so that the run starts where the abstraction does.
v0 = initial_velocity === nothing ? (x_lb[1] + x_ub[1]) / 2 : sum(initial_velocity) / 2
x0 = MP.get_coord_by_pos(
    disc.state_grid,
    MP.get_pos_by_coord(disc.state_grid, SVector(v0, q0...)),
)
println("x0 = ", x0, "  (ṗ₁ at x0: ", lgb.dynamics.f(x0, SVector(0.0, 0.0, 0.0))[1], ")")
println("swing foot at x0: ", RP.swing_foot_position(SVector{4}(x0[2:5]), geometry))
MP.get_pos_by_coord(disc.state_grid, x0) in removed && error(
    "the start cell is carved out (the swing foot at x0 is within the " *
    "$(round(100 * margin; digits = 1)) cm margin of the obstacle or below the " *
    "ground): refine the grid",
)

# All 3³ = 27 combinations of the three velocity commands by default: keeping the ankle balanced
# while the swing leg moves needs coordinated joints. Set `ONE_JOINT_PER_STEP = true` before
# including this file to allow one moving motor joint per step (7 inputs, a smaller automaton).
one_joint_per_step = @isdefined(ONE_JOINT_PER_STEP) ? ONE_JOINT_PER_STEP : false
extra_attributes =
    one_joint_per_step ?
    (; state_input_filter = (x, u) -> count(v -> abs(v) > 1e-9, u) <= 1) : (;)

# `build_abstraction` runs the Dionysos `UniformGridAbstraction` optimizer with the reach sets of
# `lgb`; the returned optimizer has the abstraction cached.
@time optimizer = LGA.build_abstraction(
    lgb;
    state_grid = disc.state_grid,
    input_grid = disc.input_grid,
    print_level = 1,
    use_implicit_mapping = true,
    execution_backend = SY.ThreadedBackend(),
    extra_attributes...,
)
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
println("transitions: ", SY.get_n_transitions(abstract_system))
println("local growth bound statistics: ", LGA.stats(lgb))

# ------------------------------------------------------------
# 4) One footstep: reach the target foothold without falling, with any velocity
# ------------------------------------------------------------

# Reach the foothold from `q0`, with any velocity at the start and at the end (by default): the
# over-approximation is nondeterministic, so the synthesis works against the worst case (no
# deterministic-only cost is given).
step_pb = RP.step_problem(
    concrete_system,
    disc.state_grid;
    q0 = q0,
    foothold = foothold,
    geometry = geometry,
    initial_velocity = initial_velocity,
    target_velocity = target_velocity,
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
@time MOI.optimize!(optimizer)
println("footstep success: ", MOI.get(optimizer, MOI.RawOptimizerAttribute("success")))

controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

# The abstraction is a certified over-approximation (up to the boundary convention of
# `LocalGrowthBound`), so a successful synthesis is guaranteed; the concrete closed loop below
# is the end-to-end check.
reached(x) = x ∈ step_pb.target_set
traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    controller,
    x0,
    400;
    stopping = reached,
)
xs = collect(ST.states(traj))
us = collect(ST.inputs(traj))
println("steps: ", length(xs) - 1, ", reached: ", reached(xs[end]))
reached(xs[end]) ||
    length(xs) > 400 ||
    println("stopped at x = ", xs[end], ": outside the controller's domain")
println(
    "ankle along the run: q₁ ∈ ",
    extrema(x[2] for x in xs),
    ", p₁ ∈ ",
    extrema(x[1] for x in xs),
)

# ------------------------------------------------------------
# 5) Visualization
# ------------------------------------------------------------

# Define `RUN_PLOTS = false` before including this file to skip it. `RP.plot_footstep` also draws
# nothing when the controller has no action at `x0` already (the default outcome of the experiment
# without a controller): `traj` has no input then.
if @isdefined(RUN_PLOTS) ? RUN_PLOTS : true
    RP.plot_footstep(
        traj,
        disc.state_grid,
        removed;
        x_lb,
        x_ub,
        tstep = disc.tstep,
        geometry,
        obstacle,
        foothold,
        file_prefix = joinpath(@__DIR__, "biped_4d_normal_form_$(scenario_name)"),
        ankle_label = "p₁ (N m s)",
        title = "Underactuated four-link biped in normal form — footstep ($(scenario_name))",
    )
end # run_plots
