# # Velocity-controlled biped: a certified footstep
#
# This example walks through the **4-D velocity-controlled biped** end to end:
# the model and its assumptions, the exact-lattice discretization, how physical
# constraints (ground, obstacles) enter the domain, the reachability objective,
# and the input slew-rate constraint that makes the synthesized command
# trackable by real motors.
#
# ## The model and its assumptions
#
# The robot is a planar biped: two legs (hip + knee joints each), the stance
# foot pinned at the origin. The **state** is the vector of the four joint
# angles `x = (θ₁, θ₂, θ₃, θ₄)` (stance hip, stance knee, swing hip, swing
# knee), and every Cartesian quantity (hip, knees, swing foot) follows by
# forward kinematics. The **input** is the vector of commanded joint
# *velocities*, and the dynamics are the pure integrator
# ```math
# \dot{x} = u .
# ```
#
# This is a deliberate change of abstraction level compared to the 6-D/8-D
# torque-controlled models of the same robot: instead of certifying the
# multibody dynamics, we **delegate the torque level to the motors' low-level
# velocity controllers** and certify the kinematic layer above them. The
# assumptions are:
#
# 1. the low-level controllers track the commanded velocity within the sampling
#    period (their tracking error is not certified here; a bounded error `ε`
#    can later be re-imported as a disturbance, degrading the abstraction from
#    exact to alternating simulation with a quantified margin);
# 2. the motion is quasi-static (slow enough that balance is not part of the
#    specification) and the stance foot stays pinned.
#
# The payoff is dramatic: the state dimension halves (angles only — velocities
# are now inputs), and `f(x, u) = x + τ·u` is exact, with no stiff multibody
# integration per transition.

using StaticArrays
using Dionysos
import LazySets
import Plots
import MathOptInterface as MOI

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction
const OPDS = DI.Optim.DiscreteSystems

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "BipedRobot",
        "4D_model",
        "robot_problem.jl",
    ),
)
import .RobotProblem as RP

geometry = RP.default_geometry()

# ## The exact-lattice discretization
#
# The discretization is chosen **commensurable**: with state-grid step `dx`,
# input-grid step `du` and sampling time `τ` such that `τ · du = dx`, every
# admissible input translates the state grid *exactly onto itself*. The
# abstraction is then deterministic and **exact** — a bisimulation, with zero
# spurious transitions — so the `CENTER_SIMULATION` approximation mode
# (unsound for general dynamics: it propagates only the cell center) becomes
# legitimate. The factory asserts the property with `MP.is_lattice_exact` and
# checks the one-cell-per-step bound `MP.intersample_safe_time_step`:

disc = RP.default_discretization(; dx = 0.1, tstep = 0.1, speed_levels = 1)
MP.is_lattice_exact(disc.state_grid, disc.input_grid, disc.tstep)

# With `speed_levels = 1` each joint velocity takes one of three values
# (`-du`, `0`, `+du`), so a step moves the state by at most one cell per axis.

# ## Physical constraints: from Cartesian space to the domain
#
# Two constraints live in the *Cartesian* space of the swing foot: the ground
# (`y ≥ 0`) and an obstacle (here a triangular step; the model also accepts a
# vector of sets — e.g. a step plus a ceiling, forcing the foot through a
# window). They are pulled back to joint space by removing every grid cell in
# which *some* configuration *might* place the swing foot inside the obstacle.
# The test is sound thanks to a Lipschitz bound of the forward kinematics: the
# foot of any configuration in a cell stays within
# `dev = Σᵢ gᵢ hᵢ/2` of the center's foot, so testing the center against the
# obstacle inflated by `dev` certifies the whole cell (`RP.carve_domain`).
#
# Together with the one-cell-per-step bound this closes the **inter-sample
# gap**: between two samples the joint-space segment stays inside the union of
# two certified cells, so the foot cannot cross the obstacle *between* samples
# either — no intermediate-time check is needed. (Bounding `τ` by the obstacle
# width instead — a natural first idea — is not sound: it prevents stepping
# fully over, but still allows grazing a corner between two free endpoints.)
#
# The margin is also an honest feasibility detector. At `dx = 0.1` it is
# `dev ≈ 5.5` cm — carving then provably *disconnects* the free space around
# an obstacle placed under the swing path: no controller exists at that
# resolution, and earlier "the trajectory jumps over the obstacle" symptoms
# came from an unsound center-only test masking exactly this. At `dx = 0.05`
# the margin halves and stepping over an 8 cm × 3 cm obstacle becomes feasible
# — that run (≈ 2.5 M cells, a few minutes) and three harder variants (the
# 16 cm × 5 cm step, a double bump, a step-plus-ceiling *window*) live in
# `examples/BipedRobot/biped_4d_velocity.jl`. Here we keep the documentation
# build light: coarse grid, and the obstacle placed *behind* the start, off
# the swing path.

obstacle =
    LazySets.VPolygon([SVector(-0.39, 0.0), SVector(-0.35, 0.02), SVector(-0.31, 0.0)])

x_bar = 0.6
domain = RP.RobotDomainConfig(;
    x_lb = SVector(-x_bar, -x_bar, -x_bar, -x_bar),
    x_ub = SVector(x_bar, x_bar, x_bar, x_bar),
    u_lb = SVector(-disc.u_max, -disc.u_max, -disc.u_max, -disc.u_max),
    u_ub = SVector(disc.u_max, disc.u_max, disc.u_max, disc.u_max),
)

concrete_system = RP.system(;
    tstep = disc.tstep,
    domain = domain,
    obstacle = obstacle,
    state_grid = disc.state_grid,
    geometry = geometry,
)
concrete_system.X isa UT.SetMinus

# ## Building the abstraction

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    PR.AlternatingSimulationProblem(concrete_system, nothing),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), disc.state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), disc.input_grid)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
(SY.get_n_state(abstract_system), SY.get_n_input(abstract_system))

# ## The objective: one footstep
#
# The task is a reachability specification: from the initial configuration
# `x0`, reach *some* configuration that places the swing foot on the target
# foothold. That target is a 2-D manifold in the 4-D joint space — the model
# enumerates it by sweeping the stance angles and solving the swing leg's
# two-link inverse kinematics (`RP.target_set`), then quantizing to grid
# cells.

x0 = SVector(0.2, 0.0, -0.2, 0.0)
foothold = SVector(0.2, 0.0)
step_pb = RP.step_problem(
    concrete_system,
    disc.state_grid;
    x0 = x0,
    foothold = foothold,
    geometry = geometry,
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
MOI.optimize!(optimizer)
MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

# Closed loop on the discrete-time system:

controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

reached(x) = x ∈ step_pb.target_set
traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    controller,
    x0,
    100;
    stopping = reached,
)
(length(collect(ST.states(traj))) - 1, reached(collect(ST.states(traj))[end]))

# ## The acceleration constraint
#
# Because the inputs are velocities, bounding the difference between
# *consecutive* inputs is an **acceleration (slew-rate) limit** — precisely
# what makes the command trackable by the motors, closing the loop with
# assumption 1. This is the classical turn-restricted shortest path: the value
# function lives on (state, input) *pairs* (the line graph), and the resulting
# controller is dynamic, with the previously played input as its memory
# (`OPDS.BoundedInputVariation`). Here:
# at most one speed notch of change per joint per step, starting from and
# ramping down to rest.

rest = SVector(0.0, 0.0, 0.0, 0.0)
slew = OPDS.BoundedInputVariation(
    (u1, u2) -> maximum(abs.(u1 - u2)),
    disc.u_max;
    target_input = rest,
    initial_input = rest,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("bounded_input_variation"), slew)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
MOI.optimize!(optimizer)
MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

#-

slew_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
slew_traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    slew_controller,
    x0,
    100;
    stopping = reached,
)
slew_us = collect(ST.inputs(slew_traj))
maximum(maximum(abs.(slew_us[k + 1] - slew_us[k])) for k in 1:(length(slew_us) - 1))

# The measured maximum input variation sits at the bound — the constraint is
# active and satisfied along the whole run.

# ## Visualizing the footstep

slew_xs = collect(ST.states(slew_traj))
plot_robot = RP.system_plot!(;
    geometry = geometry,
    obstacle = obstacle,
    foothold = foothold,
    xlims = (-0.6, 0.6),
    ylims = (-0.05, 0.5),
)
fig = Plots.plot(; aspect_ratio = :equal, legend = false)
for (k, x) in enumerate(slew_xs)
    (k % 4 == 1 || k == length(slew_xs)) && plot_robot(fig, x, rest)
end
fig

# The blue stance leg stays planted, the red swing leg travels from behind to
# the green foothold, and every visited configuration — including the segments
# *between* samples — is certified to keep the swing foot on the free side of
# the carved cells.
#
# ## What is certified, and what is not
#
# **Certified**: the closed-loop joint trajectory reaches the target foothold
# manifold, never drives the swing foot into an obstacle or below the ground
# (at sampling instants *and* in between), and respects the input slew-rate
# limit with rest-to-rest ramps.
#
# **Assumed**: the low-level velocity tracking (assumption 1) and the
# quasi-static regime (assumption 2). Quantifying the tracking error and
# re-importing it as a bounded disturbance is the natural next step — the
# growth-bound machinery accepts it directly, at the price of an abstraction
# that is no longer exact but remains sound.
