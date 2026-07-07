using StaticArrays, MathematicalSystems, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/Pendulum/double_pendulum.jl");

# concrete_problem = DoublePendulum.safety_problem(; objective = "safety_down")
concrete_problem = DoublePendulum.optimal_control_problem()

hx = SVector(5*(pi/180.0), 5*(pi/180.0), 0.5, 0.5)

u0 = SVector(0.0)
hu = SVector(0.3)

periodic_dims = SVector(1, 2)
periods = SVector(2*pi, 2*pi)
periodic_start = SVector(-pi, -pi)

Δt = 0.1

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), MP.GridFree(u0, hu))
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    SimplePendulum.jacobian_bound(),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.CENTER_SIMULATION,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> ST.NewFastIndexedAutomatonList(n, m),
)

MOI.optimize!(optimizer);

# Get the results
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")
total_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

nstep = 100
target_set =
    UT.set_in_period(concrete_problem.target_set, periodic_dims, periods, periodic_start)
nstep = 100
reached(x) = x ∈ target_set

zero_controller = MathematicalSystems.BlackBoxMap(4, 1, x -> SVector(0.0))
x0 = SVector(0.0, 0.0, 0.0, 0.0)
x_traj, u_traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller, # zero_controller
    x0,
    nstep;
    stopping = reached,
    wrap = ST.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# Here we display the coordinate projection on the two first components of the state space along the trajectory.
fig = plot(; aspect_ratio = :equal);
plot!(
    UT.set_in_period(concrete_problem.system.X, periodic_dims, periods, periodic_start);
    color = :grey,
    opacity = 1.0,
    label = "",
);
plot!(
    UT.set_in_period(concrete_problem.target_set, periodic_dims, periods, periodic_start);
    color = :red,
    opacity = 0.4,
    label = "Target set",
);
plot!(
    UT.set_in_period(concrete_problem.initial_set, periodic_dims, periods, periodic_start);
    color = :green,
    opacity = 0.8,
    label = "Initial set",
);
plot!(x_traj; ms = 2.0, arrows = false)
display(fig)

# ------------------------------------------------------------
# Animation with dashboard
# ------------------------------------------------------------

system_plot! = DoublePendulum.system_plot!()
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    x_traj,
    u_traj;
    xdims = (1, 2),
    udims = (1,),
    Δt = Δt,
    fps = 5,
    # filename = "double_pendulum_dashboard.mp4",
    xlabel_state = "θ₁ [rad]",
    ylabel_state = "θ₂ [rad]",
    xlabel_input = "time [s]",
    ylabel_input = "τ",
    xlims_state = (-π, π),
    ylims_state = (-π, π),
)

# ------------------------------------------------------------
# Visualization with RigidBodyDynamics
# ------------------------------------------------------------

# using RigidBodyDynamics
# using MeshCat
# using MeshCatMechanisms

# urdf = joinpath(
#     dirname(dirname(pathof(Dionysos))),
#     "problems",
#     "Pendulum",
#     "DoublePendulum.urdf",
# )

# Dionysos.animate_mechanism_trajectory(
#     urdf,
#     x_traj;
#     joint_names = ["joint1", "joint2"],
#     configuration = x -> [
#         x[1],          # joint1 angle
#         x[2] - x[1],   # joint2 relative angle
#     ],
#     Δt = Δt,
# )
