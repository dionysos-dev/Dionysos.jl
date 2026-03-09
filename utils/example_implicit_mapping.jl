using StaticArrays, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const SC = OP.SymbolicCertifier

include("../problems/Pendulum/simple_pendulum.jl");

concrete_problem =
    SimplePendulum.optimal_control_problem(; objective = "reachability_up_low_power") # reachability_up_high_power, reachability_up_medium_power, reachability_up_low_power, _O_ = nothing
transition_cost(x, u) = exp(10*abs(u[1]))
# concrete_problem.transition_cost = transition_cost

x0 = SVector(0.0, 0.0)
hx = SVector(3*(pi/180.0), 0.05)

u0 = SVector(0.0)
hu = SVector(0.3)

periodic_dims = SVector(1)
periods = SVector(2*pi)
periodic_start = SVector(0.5) # SVector(-pi)

tstep = 0.1

state_grid = MP.GridFree(x0, hx)
state_grid = MP.get_grid_in_periods(periodic_dims, periods, periodic_start, hx)
mapping_region = UT.HyperRectangle(SVector(-2pi, -8.0), SVector(2.5*pi, 8.0))

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

# MOI.set(optimizer, MOI.RawOptimizerAttribute("XMapping"), XMapping) or:
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("mapping_region"), mapping_region)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)

# MOI.set(optimizer, MOI.RawOptimizerAttribute("Xset"), Xset) or:
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_stateset"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("incl_mode"), MP.OUTER)

# MOI.set(optimizer, MOI.RawOptimizerAttribute("Rset"), Rset)

MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), MP.GridFree(u0, hu))
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    SimplePendulum.jacobian_bound(),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # GROWTH, CENTER_SIMULATION
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), true)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.NewIndexedAutomatonList(n, m),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("sparse_input"), true) ## 
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer);

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"));
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"));
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"));
abstract_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"));

abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")
total_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

# ------------------------------------------------------------
# Closed loop simulation
# ------------------------------------------------------------

target_set =
    UT.set_in_period(concrete_problem.target_set, periodic_dims, periods, periodic_start)
nstep = 100
reached(x) = x ∈ target_set

x0 = SVector(UT.sample(concrete_problem.initial_set)...)
x_traj, u_traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
    wrap = ST.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
Rset = MOI.get(optimizer, MOI.RawOptimizerAttribute("Rset"))
Xset = MOI.get(optimizer, MOI.RawOptimizerAttribute("Xset"))
XMapping = MOI.get(optimizer, MOI.RawOptimizerAttribute("XMapping"))

fig = plot(; aspect_ratio = :equal);
concrete_system = concrete_problem.system
plot!(XMapping; color = :blue, opacity = 0.15, label = "Rset")
# plot!((Rset, XMapping); color = :grey, opacity = 0.15, label = "Rset")
# plot!((Xset, XMapping); color = :grey, opacity = 0.25, label = "Xset")
plot!(
    UT.set_in_period(concrete_system.X, periodic_dims, periods, periodic_start);
    color = :grey,
    hole_color = :black,
    opacity = 1.0,
    label = "",
);
plot!(
    UT.set_in_period(concrete_problem.initial_set, periodic_dims, periods, periodic_start);
    color = :green,
    opacity = 0.2,
    label = "Initial set",
);
plot!(target_set; color = :red, opacity = 0.8, label = "Target set");
plot!(x_traj; ms = 2.0, arrows = false)
display(fig)
