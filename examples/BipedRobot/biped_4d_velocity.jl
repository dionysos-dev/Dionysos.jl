# # 4-D velocity-controlled biped: one certified footstep over an obstacle
#
# The biped is controlled at the *velocity* level: states are the 4 joint
# angles, inputs the joint angular velocities (`ẋ = u`), and the motors' low
# level controllers track the commanded velocities. Compared to the 6D/8D
# torque models this halves the abstraction dimension and makes the dynamics a
# translation.
#
# Three ingredients make the run below fully certified:
#
# 1. **Exact lattice**: `tstep · du = dx`, so every input translates the state
#    grid onto itself — the abstraction is deterministic and *exact* (a
#    bisimulation), which is what legitimizes `CENTER_SIMULATION`
#    (`MP.is_lattice_exact` is asserted by the model's discretization factory).
# 2. **Sound obstacle carving**: the Cartesian obstacle (a step in front of the
#    stance foot that the swing foot must climb over) is pulled back to joint
#    space with a Lipschitz margin — a kept cell provably keeps the swing foot
#    clear (`RobotProblem.carve_domain`).
# 3. **Swept-cell transitions** (`MP.cells_on_segment`): the inter-sample
#    trajectory of the integrator is the straight joint-space segment, and a
#    transition is kept only when *every* grid cell that segment crosses is
#    admissible — so the swing foot cannot cross the obstacle *between*
#    samples either. This lifts the one-cell-per-step speed cap: the input
#    alphabet has two speed levels per joint (`±0.5, ±1.0` rad/s), and the
#    fast inputs move two cells per step.
#
# Resolution matters: at `dx = 0.1` the Lipschitz margin (~5.5 cm) provably
# disconnects the free space — no controller exists. At `dx = 0.05` the margin
# drops to ~2.7 cm and the step becomes feasible. The run takes a few minutes
# (≈ 2.5 M cells); inputs are restricted to one joint per step to keep the
# automaton laptop-sized.
#
# The second synthesis pass adds the **input slew-rate constraint**
# (`OPDS.BoundedInputVariation`): consecutive joint-velocity commands differ by
# at most one speed notch and the motion starts from and ramps down to rest —
# the profile a real motor controller can track.

using StaticArrays
using MathematicalSystems
using Dionysos
import MathOptInterface as MOI
import LazySets
import Plots

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const OPDS = OP.DiscreteSystems

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

# ------------------------------------------------------------
# 1) Scenario
# ------------------------------------------------------------
# Pick with the BIPED_SCENARIO environment variable. All four are BFS-verified
# connected under sound carving at dx = 0.05:
#   step        — 8 cm × 3 cm step, |θ| ≤ 1.2 (the default; ~52 steps)
#   riccardo    — the original 16 cm × 5 cm step: infeasible on a symmetric
#                 domain, feasible once the swing leg gets ±1.4 rad (~58 steps)
#   wall        — a thin 10 cm wall: the certified clearance (~12.7 cm with
#                 the Lipschitz margin) is close to the kinematic maximum
#   window      — step + ceiling at 10 cm: the foot threads a certified
#                 vertical slot (~52 steps)

_triangle(cx, hw, h) =
    LazySets.VPolygon([SVector(cx - hw, 0.0), SVector(cx, h), SVector(cx + hw, 0.0)])

const SCENARIOS = Dict(
    :step => (;
        obstacles = [_triangle(0.0, 0.04, 0.03)],
        x_lb = SVector(-1.2, -1.2, -1.2, -1.2),
        x_ub = SVector(1.2, 1.2, 1.2, 1.2),
        foothold = SVector(0.2, 0.0),
    ),
    :riccardo => (;
        obstacles = [_triangle(0.0, 0.08, 0.05)],
        x_lb = SVector(-0.8, -0.8, -1.4, -1.4),
        x_ub = SVector(0.8, 0.8, 1.4, 1.4),
        foothold = SVector(0.2, 0.0),
    ),
    :wall => (;
        obstacles = [_triangle(0.0, 0.02, 0.10)],
        x_lb = SVector(-1.2, -1.2, -1.2, -1.2),
        x_ub = SVector(1.2, 1.2, 1.2, 1.2),
        foothold = SVector(0.2, 0.0),
    ),
    :window => (;
        obstacles = [
            _triangle(0.0, 0.04, 0.03),
            LazySets.Hyperrectangle(;
                low = SVector(-0.05, 0.10),
                high = SVector(0.05, 0.30),
            ),
        ],
        x_lb = SVector(-1.2, -1.2, -1.2, -1.2),
        x_ub = SVector(1.2, 1.2, 1.2, 1.2),
        foothold = SVector(0.2, 0.0),
    ),
)

scenario_name = Symbol(get(ENV, "BIPED_SCENARIO", "step"))
scenario = SCENARIOS[scenario_name]
println("scenario: ", scenario_name)

geometry = RP.default_geometry()

# Two speed levels per joint (`u ∈ {-1, -0.5, 0, 0.5, 1}` rad/s): steps of up
# to two cells per axis, made sound by the swept-cell transition validation
# below (`MP.swept_input_filter`).
disc = RP.default_discretization(;
    dx = 0.05,
    tstep = 0.1,
    speed_levels = 2,
    swept_transitions = true,
)

obstacle = scenario.obstacles
domain = RP.RobotDomainConfig(;
    x_lb = scenario.x_lb,
    x_ub = scenario.x_ub,
    u_lb = SVector(-disc.u_max, -disc.u_max, -disc.u_max, -disc.u_max),
    u_ub = SVector(disc.u_max, disc.u_max, disc.u_max, disc.u_max),
)

x0 = SVector(0.2, 0.0, -0.2, 0.0)
foothold = scenario.foothold

println("Carving the obstacle out of the joint-space domain…")
X_box = LazySets.Hyperrectangle(; low = scenario.x_lb, high = scenario.x_ub)
removed = RP.infeasible_cells(X_box, disc.state_grid, obstacle, geometry)
concrete_system = RP.system(;
    tstep = disc.tstep,
    domain = domain,
    state_grid = disc.state_grid,
    removed_cells = removed,
)

# ------------------------------------------------------------
# 2) Abstraction (exact lattice, one joint per step)
# ------------------------------------------------------------

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
# Two combined restrictions on (state, input) pairs:
# - one joint per step, which keeps the automaton at 17 effective inputs
#   (≈ 40 M transitions) instead of 5⁴ = 625 (beyond laptop memory);
# - the swept-cell check: a multi-cell step is kept only when every grid cell
#   crossed by its inter-sample segment is admissible — the sound replacement
#   for the one-cell-per-step speed cap.
swept = MP.swept_input_filter(disc.state_grid, disc.tstep, removed)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("state_input_filter"),
    (x, u) -> count(v -> abs(v) > 1e-9, u) <= 1 && swept(x, u),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("intersample_checked"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("execution_backend"), SY.ThreadedBackend())
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

@time MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))

# ------------------------------------------------------------
# 3) One footstep: reach the target foothold
# ------------------------------------------------------------

# The (x, u) -> 1 cost routes the synthesis to the deterministic Dijkstra,
# which needs no state × input counter matrix (the uniform-cost BFS would
# allocate one over all 3⁴ input symbols).
step_pb = RP.step_problem(
    concrete_system,
    disc.state_grid;
    x0 = x0,
    foothold = foothold,
    geometry = geometry,
)
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

reached(x) = x ∈ step_pb.target_set
traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    controller,
    x0,
    400;
    stopping = reached,
)
xs = collect(ST.states(traj))
println("steps: ", length(xs) - 1, ", reached: ", reached(xs[end]))

# ------------------------------------------------------------
# 4) Same footstep with the acceleration (slew-rate) limit
# ------------------------------------------------------------

# Consecutive velocity commands within one speed notch (`du`, not `u_max`!)
# per joint, starting from and ramping down to rest: reaching full speed takes
# two steps (0 → 0.5 → 1.0 rad/s) and reversals must ramp back through zero.
rest = SVector(0.0, 0.0, 0.0, 0.0)
slew = OPDS.BoundedInputVariation(
    (u1, u2) -> maximum(abs.(u1 - u2)),
    disc.du;
    target_input = rest,
    initial_input = rest,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("bounded_input_variation"), slew)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), step_pb)
@time MOI.optimize!(optimizer)
println(
    "footstep success (slew-rate limited): ",
    MOI.get(optimizer, MOI.RawOptimizerAttribute("success")),
)

slew_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
slew_traj = ST.get_closed_loop_trajectory(
    discrete_time_system,
    slew_controller,
    x0,
    400;
    stopping = reached,
)
slew_xs = collect(ST.states(slew_traj))
slew_us = collect(ST.inputs(slew_traj))
println(
    "steps (slew-rate limited): ",
    length(slew_xs) - 1,
    ", reached: ",
    reached(slew_xs[end]),
)
max_slew =
    maximum(maximum(abs.(slew_us[k + 1] - slew_us[k])) for k in 1:(length(slew_us) - 1))
println("max input variation along the run: ", max_slew, " (bound: ", disc.du, ")")
println("max input magnitude along the run: ", maximum(maximum(abs.(u)) for u in slew_us))

# ------------------------------------------------------------
# 5) Visualization
# ------------------------------------------------------------

plot_robot = RP.system_plot!(;
    geometry = geometry,
    obstacle = obstacle,
    foothold = foothold,
    xlims = (-0.6, 0.6),
    ylims = (-0.05, 0.55),
)

png_file = joinpath(@__DIR__, "biped_4d_footstep_$(scenario_name).png")
gif_file = joinpath(@__DIR__, "biped_4d_footstep_$(scenario_name).gif")

fig = Plots.plot(; aspect_ratio = :equal, legend = false)
for (k, x) in enumerate(slew_xs)
    k % 5 == 1 || k == length(slew_xs) || continue
    plot_robot(fig, x, k <= length(slew_us) ? slew_us[k] : rest)
end
Plots.savefig(fig, png_file)
println("saved ", png_file)

# State panel = the swing angles (θ3, θ4) — where the obstacle preimage lives.
# The carved region is 4-D, so drawing its raw projection would be misleading;
# instead each frame shows its *slice*: the (θ3, θ4) cells that are blocked
# given the current stance angles (θ1, θ2).
h3, h4 = MP.get_h(disc.state_grid)[3], MP.get_h(disc.state_grid)[4]
k3_range = round(Int, scenario.x_lb[3] / h3):round(Int, scenario.x_ub[3] / h3)
k4_range = round(Int, scenario.x_lb[4] / h4):round(Int, scenario.x_ub[4] / h4)
blocked_slice! = function (p_state, x)
    pos = MP.get_pos_by_coord(disc.state_grid, x)
    shapes = Plots.Shape[]
    for k3 in k3_range, k4 in k4_range
        (pos[1], pos[2], k3, k4) in removed || continue
        c = MP.get_coord_by_pos(disc.state_grid, (pos[1], pos[2], k3, k4))
        push!(
            shapes,
            Plots.Shape(
                [c[3] - h3 / 2, c[3] + h3 / 2, c[3] + h3 / 2, c[3] - h3 / 2],
                [c[4] - h4 / 2, c[4] - h4 / 2, c[4] + h4 / 2, c[4] + h4 / 2],
            ),
        )
    end
    isempty(shapes) ||
        Plots.plot!(p_state, shapes; color = :black, alpha = 0.35, lw = 0, label = "")
    return p_state
end

DI.animate_trajectory_dashboard(
    plot_robot,
    slew_traj;
    Δt = disc.tstep,
    xdims = (3, 4),
    udims = (3, 4),
    xlabel_state = "θ3 (swing hip)",
    ylabel_state = "θ4 (swing knee)",
    xlabel_input = "u3 (swing hip)",
    ylabel_input = "u4 (swing knee)",
    xlims_state = (scenario.x_lb[3], scenario.x_ub[3]),
    ylims_state = (scenario.x_lb[4], scenario.x_ub[4]),
    state_background! = blocked_slice!,
    title = "4-D velocity-controlled biped — certified footstep ($(scenario_name))",
    filename = gif_file,
)
println("saved ", gif_file)
