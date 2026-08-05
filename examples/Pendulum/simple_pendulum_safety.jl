using StaticArrays, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "Pendulum",
        "simple_pendulum.jl",
    ),
);

concrete_problem = SimplePendulum.safety_problem(; objective = "safety_up") # safety_up, safety_down

hx = SVector(3*(pi/180.0), 0.05)

u0 = SVector(0.0)
hu = SVector(0.3)

periodic_dims = SVector(1)
periods = SVector(2*pi)
periodic_start = SVector(-pi)

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
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.FastIndexedAutomatonList(n, m),
)

MOI.optimize!(optimizer);

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"));
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"));
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"));
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"));
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")
total_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

nstep = 300

x0 = SVector(UT.sample(concrete_problem.initial_set)...)
traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
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
    UT.set_in_period(concrete_problem.safe_set, periodic_dims, periods, periodic_start);
    color = :red,
    opacity = 0.4,
    label = "Safe set",
);
plot!(
    UT.set_in_period(concrete_problem.initial_set, periodic_dims, periods, periodic_start);
    color = :green,
    opacity = 0.8,
    label = "Initial set",
);
plot!(traj; ms = 2.0, arrows = false)
display(fig)

# ------------------------------------------------------------
# Animation with dashboard
# ------------------------------------------------------------

system_plot! = SimplePendulum.system_plot!()
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1, 2),      # phase plot θ vs ω
    udims = (1,),        # input over time
    Δt = Δt,
    fps = 5,
    # filename = "simple_pendulum_dashboard.mp4",
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    xlabel_input = "time [s]",
    ylabel_input = "τ [Nm]",
    xlims_state = (-π, π),
    ylims_state = (-8, 8),
)
