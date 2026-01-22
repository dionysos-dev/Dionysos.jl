using StaticArrays, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include("../problems/simple_pendulum.jl");

concrete_problem = SimplePendulum.optimal_control_problem()

hx = SVector(pi/100.0, pi/100.0)

u0 = SVector(0.0)
hu = SVector(0.3)

periodic_dims = SVector(1)
periods = SVector(2*pi)
periodic_start = SVector(-pi)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("input_grid"),
    Dionysos.Domain.GridFree(u0, hu),
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    SimplePendulum.jacobian_bound(),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_domain"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), true) # true
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
# MOI.set(optimizer, MOI.Silent(), true)

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

target_set = concrete_problem.target_set
target_set_in_periodic =
    UT.set_in_period(target_set, periodic_dims, periods, periodic_start)

nstep = 100
function reached(x)
    if x ∈ target_set_in_periodic
        return true
    else
        return false
    end
end

x0 = SVector(UT.sample(concrete_problem.initial_set)...)
x_traj, u_traj = Dionysos.System.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
    wrap = ST.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# Here we display the coordinate projection on the two first components of the state space along the trajectory.
fig = plot(; aspect_ratio = :equal, legend = false);
# We display the concrete domain
state_space = concrete_system.X
system_domain_in_periodic =
    UT.set_in_period(state_space, periodic_dims, periods, periodic_start)
plot!(system_domain_in_periodic; color = :grey, opacity = 1.0, label = "");

# We display the abstract domain with worst-case cost
# plot!(abstract_system; value_function = abstract_value_function);

# We display the concrete specifications
plot!(concrete_problem.initial_set; color = :green, opacity = 0.2, label = "Initial set");
# plot!(target_set; dims = [1, 2], color = :red, opacity = 0.5, label = "Target set");
plot!(
    target_set_in_periodic;
    dims = [1, 2],
    color = :red,
    opacity = 0.8,
    label = "Target set in periodic domain",
);

# We display the abstract specifications
# plot!(
#     Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.initial_set);
#     color = :green,
# );
# plot!(
#     Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.target_set);
#     color = :red,
# );

# We display the concrete trajectory
plot!(x_traj; ms = 2.0, arrows = false)
