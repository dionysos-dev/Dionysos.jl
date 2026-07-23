using StaticArrays, Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/PathPlanning/path_planning.jl")

####################################################################################
########### PART 1 : Construct the classical abstraction using any mode: ###########
##### USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION  ##########
####################################################################################
transition_cost(x, u) = 1.0 + 1 / (1e-4 + (x[3] - π/2)^6) + 1 / (1e-4 + (x[3] + π/2)^6)
# penalize x[3] from being close to π/2 and -π/2, i.e., going up and down

concrete_problem = PathPlanning.problem(; simple = true, transition_cost = transition_cost)
concrete_system = concrete_problem.system

x0 = SVector(0.0, 0.0, 0.0)
h = SVector(0.2, 0.2, 0.2)
state_grid = MP.GridFree(x0, h)

u0 = SVector(0.0, 0.0)
h = SVector(0.3, 0.3)
input_grid = MP.GridFree(u0, h)

using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    PathPlanning.jacobian_bound(),
)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # USER_DEFINED GROWTH LINEARIZED CENTER_SIMULATION RANDOM_SIMULATION
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.IndexedAutomatonList(n, m),
)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"))
concrete_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_value_function"))

nstep = 300
reached(x) = x ∈ concrete_problem.target_set

x0 = SVector(0.4, 0.4, 0.0)
println(
    "Worst-case (upper bound) value for the initial point: ",
    concrete_value_function(x0),
)
traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

Xmap = SY.get_state_mapping(abstract_system)
fig = plot(; aspect_ratio = :equal)
plot!(concrete_system.X; color = :grey, opacity = 1.0, label = "")
plot!(abstract_system; value_function = abstract_value_function)
plot!(
    UT.project_set(concrete_problem.initial_set, [1, 2]);
    color = :green,
    opacity = 0.2,
    label = "Initial set",
)
plot!(
    UT.project_set(concrete_problem.target_set, [1, 2]);
    color = :red,
    opacity = 0.5,
    label = "Target set",
)
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.initial_set), Xmap);
    color = :green,
    efficient = false,
)
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.target_set), Xmap);
    color = :red,
    efficient = false,
)
plot!(traj; ms = 2.0, arrows = false, color = :blue)
display(fig)

####################################################################################
### PART 2 : Construct the determinized abstraction and solve a concrete problem ###
####################################################################################

println(
    "Is the original abstract system deterministic ? ",
    SY.is_deterministic(abstract_system),
)
determinized_abstract_system = SY.determinize_symbolic_model(abstract_system)
println(
    "Is the determinized abstract system deterministic ? ",
    SY.is_deterministic(determinized_abstract_system),
)

println(SY.get_n_transitions(abstract_system))
println(SY.get_n_transitions(determinized_abstract_system))

concrete_problem = PathPlanning.problem(; simple = true, transition_cost = transition_cost)

using JuMP
new_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(
    new_optimizer,
    MOI.RawOptimizerAttribute("abstract_system"),
    determinized_abstract_system,
)
MOI.set(
    new_optimizer,
    MOI.RawOptimizerAttribute("discrete_time_system"),
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("sparse_input"), true)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.IndexedAutomatonList(n, m),
)

MOI.optimize!(new_optimizer)

abstract_problem = MOI.get(new_optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
concrete_controller =
    MOI.get(new_optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_value_function =
    MOI.get(new_optimizer, MOI.RawOptimizerAttribute("abstract_value_function"))
concrete_value_function =
    MOI.get(new_optimizer, MOI.RawOptimizerAttribute("concrete_value_function"))

println(
    "Heuristic (lower bound) value for the initial point: ",
    concrete_value_function(x0),
)

new_traj = ST.get_closed_loop_trajectory(
    MOI.get(new_optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

Xmap = SY.get_state_mapping(determinized_abstract_system)
fig = plot(; aspect_ratio = :equal)
plot!(concrete_system.X; color = :grey, opacity = 1.0, label = "")
plot!(determinized_abstract_system; value_function = abstract_value_function)
plot!(
    UT.project_set(concrete_problem.initial_set, [1, 2]);
    color = :green,
    opacity = 0.2,
    label = "Initial set",
)
plot!(
    UT.project_set(concrete_problem.target_set, [1, 2]);
    color = :red,
    opacity = 0.5,
    label = "Target set",
)
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.initial_set), Xmap);
    color = :green,
    efficient = false,
)
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.target_set), Xmap);
    color = :red,
    efficient = false,
)
plot!(traj; ms = 2.0, arrows = false, color = :blue)
plot!(new_traj; ms = 2.0, arrows = false, color = :red)
display(fig)
