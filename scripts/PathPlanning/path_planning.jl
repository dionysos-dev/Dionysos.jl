using StaticArrays, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/PathPlanning/path_planning.jl");

concrete_problem = PathPlanning.problem(; simple = true)
concrete_system = concrete_problem.system

x0 = SVector(0.0, 0.1, 0.0);
hx = SVector(0.2, 0.2, 0.2);
u0 = SVector(0.0, 0.0);
hu = SVector(0.3, 0.3);
Δt = 0.3;
periodic_dims = SVector(2); # SVector(1, 2);, 
periods = SVector(10.0); # SVector(4.0, 10.0);
periodic_start = SVector(0.0); # SVector(0.0, 0.0);
mapping_region = UT.box(SVector(0.0, 0.0, -pi - 0.4), SVector(4.0, 11.0, pi + 0.4))

# Intantiate the optimizer
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

# Set the control problem
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

# State Mapping attributes
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("mapping_region"), mapping_region)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)

# State set attributes
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_stateset"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("incl_mode"), MP.INNER)

# Input Mapping attributes
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), MP.GridFree(u0, hu))
# Time Mapping attributes
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)

# Other attributes
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    PathPlanning.jacobian_bound(),
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false) # true
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer);

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"));
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"));
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"));
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"));
abstract_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"));
concrete_system = concrete_problem.system
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")
total_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

## Export in csv file the controller, and reload it
# filename = "concrete_controller"
# AB.UniformGridAbstraction.export_controller_csv(optimizer, filename)
# AB.UniformGridAbstraction.import_controller_csv(filename)

target_set = concrete_problem.target_set
target_set_in_periodic =
    UT.set_in_period(target_set, periodic_dims, periods, periodic_start)

nstep = 100
reached(x) = x ∈ target_set_in_periodic

x0 = SVector(0.4, 0.4, 0.0)
x_traj, u_traj = Dionysos.System.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
    wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# Here we display the coordinate projection on the two first components of the state space along the trajectory.
fig = plot(; aspect_ratio = :equal, legend = false);
# We display the concrete domain
state_space = concrete_system.X
system_domain_in_periodic =
    UT.set_in_period(state_space, periodic_dims, periods, periodic_start)
plot!(system_domain_in_periodic; color = :grey, opacity = 1.0, label = "");

# We display the abstract domain with worst-case cost
XMapping = SY.get_state_mapping(abstract_system)
Xset = SY.get_state_set(abstract_system)
plot!(abstract_system; efficient = false, value_function = abstract_value_function);
# plot!(XMapping; value_function=abstract_value_function, efficient=false)
# plot!((Xset, XMapping); efficient=true, color = :yellow)

# We display the concrete specifications
plot!(
    UT.project_set(concrete_problem.initial_set, [1, 2]);
    color = :green,
    opacity = 0.2,
    label = "Initial set",
);
plot!(UT.project_set(target_set, [1, 2]); color = :red, opacity = 0.5, label = "Target set");
plot!(
    target_set_in_periodic;
    dims = [1, 2],
    color = :red,
    opacity = 0.8,
    label = "Target set in periodic domain",
);

# We display the abstract specifications
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.initial_set), XMapping);
    color = :green,
);
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.target_set), XMapping);
    color = :red,
    efficient = false,
);

# We display the concrete trajectory
plot!(x_traj; ms = 2.0, arrows = false)

# ------------------------------------------------------------
# Animation with dashboard
# ------------------------------------------------------------

_X_ = UT.box(SVector(0.0, 0.0, -pi - 0.4), SVector(4.0, 10.0, pi + 0.4))
obstacles = PathPlanning.get_obstacles(_X_)
system_plot! = PathPlanning.system_plot!(;
    obstacles = obstacles,
    xlims = (-0.1, 10.1),
    ylims = (-0.1, 10.1),
)
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    x_traj,
    u_traj;
    xdims = (1, 2),
    udims = (1, 2),
    Δt = Δt,
    fps = 5,
    # filename = "path_planning_dashboard.mp4",
    xlabel_state = "x",
    ylabel_state = "y",
    xlabel_input = "v",
    ylabel_input = "δ",
    xlims_state = (-0.1, 10.1),
    ylims_state = (-0.1, 10.1),
)
