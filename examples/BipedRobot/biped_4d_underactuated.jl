# # Underactuated four-link biped: one footstep over an obstacle
#
# Velocity-controlled design model of the four-link biped of
# `problems/BipedRobot/4D_model_underactuated/robot_problem.jl` (a symbolic `Multibody2D` model, the
# same robot as `four_link_biped.urdf`). Only the three
# motor joints (grounded knee `q₂`, hip `q₃`, swing knee `q₄`) are commanded, in
# velocity; the grounded ankle `q₁` is **unactuated** and obeys the first
# Euler–Lagrange row of the multibody model. The state and input are
#
#     x = [q̇₁, q₁, q₂, q₃, q₄],   u = [v₂, v₃, v₄],
#
# with `ẋ = [q̈₁, q̇₁, v₂, v₃, v₄]` (commands held constant over a sampling period, i.e. zero
# actuated accelerations) and `q̈₁ = h₁ / D₁₁` from the first Euler–Lagrange row.
#
# Compared to the fully actuated 4-D example (`biped_4d_velocity.jl`, `ẋ = u`):
#
# - the ankle is an inverted pendulum (time constant ≈ 0.2 s), so the controller
#   has to keep the robot from falling while it moves the swing foot over the
#   obstacle, using only the posture of the motor joints;
# - the dynamics are nonlinear in `(q̇₁, q₁)`: `tstep · du = dx` still keeps the
#   motor coordinates on the lattice, but the ankle coordinates are rounded to the
#   grid, so the `CENTER_SIMULATION` abstraction is **not certified** here. The
#   closed-loop trajectory of the concrete (RK4-integrated) system below is what
#   tells whether the footstep succeeds;
# - angles are **counter-clockwise**: a positive `q₁` leans the robot backwards, a positive `q₃`
#   swings the swing foot forwards and the knee bends backwards for a negative `q₄`;
# - the multi-cell ankle moves are not swept-cell checked, so the solver's
#   inter-sample warning may fire.
#
# The obstacle carving of the swing foot is still sound (Lipschitz pullback of the
# forward kinematics, compiled from the same symbolic model as the dynamics). RigidBodyDynamics is not
# needed here (it is only used by `BipedFourLink.jl`, to simulate the whole robot with torques); the
# symbolic model reproduces its dynamics to machine precision (`tests/test_biped_four_link.jl`).
# The model is the one of the normal-form experiment (`biped_4d_normal_form.jl`); only the state
# (`q̇₁` instead of `p₁`) and the abstraction differ.
#
# Run it with (the `BipedRobot` environment provides Plots and the other packages, to which
# Symbolics must be added):
#
# ```
# julia -t auto --project=problems/BipedRobot examples/BipedRobot/biped_4d_underactuated.jl
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

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "BipedRobot",
        "4D_model_underactuated",
        "robot_problem.jl",
    ),
)
import .UnderactuatedRobotProblem as RP

# ------------------------------------------------------------
# 1) Scenario: the `step` of the fully actuated example
# ------------------------------------------------------------
# An 8 cm × 3 cm step at the stance foot that the swing foot must clear, and a
# foothold 20 cm in front of the stance foot.

_triangle(cx, hw, h) =
    LazySets.VPolygon([SVector(cx - hw, 0.0), SVector(cx, h), SVector(cx + hw, 0.0)])

scenario_name = :step
obstacle = @isdefined(OBSTACLE) ? OBSTACLE : [_triangle(0.0, 0.04, 0.03)]
foothold = SVector(0.2, 0.0)

model = RP.symbolic_model()
geometry = RP.Geometry(model)                         # kinematics, from the symbolic model
dynamics = RP.VelocityControlledDynamics(; model)     # q̈₁ = h₁ / D₁₁, compiled

# Set `DISC` before including this file to try another resolution. The motor
# joints move by exactly one cell (`dx`) per step at speed `du = dx / tstep`; the
# ankle angle and velocity have their own steps `dq1`, `dq̇1`.
disc =
    @isdefined(DISC) ? DISC :
    RP.default_discretization(; dx = 0.1, dq1 = 0.05, dq̇1 = 0.1, tstep = 0.1)

# State box `[q̇₁, q₁, q₂, q₃, q₄]`: the robot "falls" when it leaves the ankle
# bounds, which the reach specification forbids. Every posture placing the swing
# foot on the foothold within these bounds has the ankle falling forwards
# (q̈₁ ≲ -4 rad/s² at rest), so the step is a controlled fall into the foothold and
# the q̇₁ range must leave room for it. The swing hip only needs to swing forwards
# (q₃ > 0) and the swing knee to bend backwards (q₄ ≤ 0, the foot behind the knee).
# Set `X_BOUNDS = (x_lb, x_ub)` before including this file to change them.
x_lb, x_ub =
    @isdefined(X_BOUNDS) ? X_BOUNDS :
    (SVector(-1.5, -0.4, -0.8, -0.4, -1.2), SVector(1.5, 0.4, 0.8, 1.2, 0.2))
domain = RP.RobotDomainConfig(;
    x_lb = x_lb,
    x_ub = x_ub,
    u_lb = SVector(-disc.u_max, -disc.u_max, -disc.u_max),
    u_ub = SVector(disc.u_max, disc.u_max, disc.u_max),
)

# The footstep starts from the configuration `q0 = [q₁, q₂, q₃, q₄]` (stance knee straight, the
# swing foot 15 cm behind the stance foot: swing hip 0.2 rad behind, swing knee bent backwards
# by 0.5 rad, which keeps the foot outside the carved neighbourhood of the obstacle) and ends on the
# foothold. **No velocity is imposed**: the controller must work from every ankle velocity of
# `INITIAL_VELOCITY` and may end with any ankle velocity of `TARGET_VELOCITY`; both are
# `(vlo, vhi)` ranges on the ankle velocity `q̇₁` and default to its whole interval of definition
# (`nothing`). The swing foot has no ground reaction force, so requiring rest is needlessly
# conservative. Set `Q0`, `INITIAL_VELOCITY` or `TARGET_VELOCITY` before including this file to
# change them, e.g. `INITIAL_VELOCITY = (-0.3, 0.3)`.
q0 = @isdefined(Q0) ? Q0 : SVector(-0.05, 0.0, -0.2, -0.5)
initial_velocity = @isdefined(INITIAL_VELOCITY) ? INITIAL_VELOCITY : nothing
target_velocity = @isdefined(TARGET_VELOCITY) ? TARGET_VELOCITY : nothing

# Concrete start of the closed-loop run below: the configuration `q0` with the ankle velocity in the
# middle of the initial range, snapped to the grid so that the run starts where the abstraction does.
v0 = initial_velocity === nothing ? (x_lb[1] + x_ub[1]) / 2 : sum(initial_velocity) / 2
x0 = MP.get_coord_by_pos(
    disc.state_grid,
    MP.get_pos_by_coord(disc.state_grid, SVector(v0, q0...)),
)
println("x0 = ", x0, "  (q̈₁ at x0: ", dynamics(x0, SVector(0.0, 0.0, 0.0))[1], ")")

println("Carving the obstacle out of the state domain…")
X_box = LazySets.Hyperrectangle(; low = x_lb, high = x_ub)
removed = RP.infeasible_cells(X_box, disc.state_grid, obstacle, geometry)

margin = RP.swing_foot_deviation_bound(geometry, MP.get_h(disc.state_grid) ./ 2)
println("Lipschitz carving margin: ", round(100 * margin; digits = 1), " cm")
MP.get_pos_by_coord(disc.state_grid, x0) in removed && error(
    "the start cell is carved out (the swing foot at x0 is within the " *
    "$(round(100 * margin; digits = 1)) cm margin of the obstacle or below the " *
    "ground): refine the grid",
)

concrete_system = RP.system(; dynamics = dynamics, domain = domain, removed_cells = removed)

# ------------------------------------------------------------
# 2) Abstraction (center simulation of the RK4 roll-out)
# ------------------------------------------------------------

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    PR.AlternatingSimulationProblem(concrete_system, nothing),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), disc.state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), disc.input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), disc.tstep)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_mapping"), true)
# All 3³ = 27 combinations of the three velocity commands by default: keeping the
# ankle balanced while the swing leg moves needs coordinated joints. Set
# `ONE_JOINT_PER_STEP = true` before including this file to allow one moving motor
# joint per step (7 inputs, a smaller automaton).
one_joint_per_step = @isdefined(ONE_JOINT_PER_STEP) ? ONE_JOINT_PER_STEP : false
if one_joint_per_step
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("state_input_filter"),
        (x, u) -> count(v -> abs(v) > 1e-9, u) <= 1,
    )
end
MOI.set(optimizer, MOI.RawOptimizerAttribute("execution_backend"), SY.ThreadedBackend())
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

@time MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))

# ------------------------------------------------------------
# 3) One footstep: reach the target foothold without falling, with any velocity
# ------------------------------------------------------------

step_pb = RP.step_problem(
    concrete_system,
    disc.state_grid;
    q0 = q0,
    foothold = foothold,
    geometry = geometry,
    initial_velocity = initial_velocity,
    target_velocity = target_velocity,
)
# The (x, u) -> 1 cost routes the synthesis to the deterministic Dijkstra.
step_pb = PR.OptimalControlProblem(
    concrete_system,
    step_pb.initial_set,
    step_pb.target_set,
    nothing,
    (x, u) -> 1.0,
    PR.Infinity(),
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
@time MOI.optimize!(optimizer)
println("footstep success: ", MOI.get(optimizer, MOI.RawOptimizerAttribute("success")))

controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

# The abstraction is not certified: this concrete closed loop is the check.
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
# The closed loop only stops early when the controller has no action for the
# current state: the concrete ankle drifted away from the rounded abstract one.
reached(xs[end]) ||
    length(xs) > 400 ||
    println("stopped at x = ", xs[end], ": outside the controller's domain")
println(
    "ankle along the run: q₁ ∈ ",
    extrema(x[2] for x in xs),
    ", q̇₁ ∈ ",
    extrema(x[1] for x in xs),
)

# ------------------------------------------------------------
# 4) Visualization
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
        file_prefix = joinpath(@__DIR__, "biped_4d_underactuated_$(scenario_name)"),
        ankle_label = "q̇₁ (rad/s)",
        title = "Underactuated four-link biped — footstep ($(scenario_name))",
    )
end # run_plots
