module TestMain
using Test     #src

using StaticArrays, Plots

# At this point, we import Dionysos.
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

# And the file defining the hybrid system for this problem
include("../../../problems/path_planning.jl")

# ### Definition of the problem

# Now we instantiate the problem using the function provided by [PathPlanning.jl](@__REPO_ROOT_URL__/problems/PathPlanning.jl) 
concrete_problem = PathPlanning.problem(; simple = true)
concrete_system = concrete_problem.system

# ### Definition of the abstraction

# Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):
x0 = SVector(0.0, 0.0, 0.0)
h = SVector(0.2, 0.2, 0.2)
state_grid = DO.GridFree(x0, h)

# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(0.0, 0.0)
h = SVector(0.3, 0.3)
input_grid = DO.GridFree(u0, h)

# We now solve the optimal control problem with the `Abstraction.UniformGridAbstraction.Optimizer`.

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
MOI.optimize!(optimizer)

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))


nstep = 150
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end
x0 = SVector(0.4, 0.4, 0.0)
control_trajectory = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

fig = plot(; aspect_ratio = :equal)
plot!(concrete_system.X; color = :yellow, opacity = 0.5, label = "")
plot!(abstract_system.Xdom; color = :blue, opacity = 0.5, label = "", efficient=false)
plot!(concrete_problem.initial_set; color = :green, label = "", opacity = 0.2)
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, label = "", opacity = 0.2)
plot!(
        SY.get_domain_from_states(abstract_system, abstract_problem.initial_set);
        color = :green,
        label = "",
    )
plot!(
        SY.get_domain_from_states(abstract_system, abstract_problem.target_set);
        color = :red,
        label = "",
    )
plot!(control_trajectory; ms = 0.5)
display(fig)



###### New 
using JuMP
new_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("abstract_system"), abstract_system)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.optimize!(new_optimizer)


println(SY.is_deterministic(abstract_system))
deterministic_abstract_system = SY.determinize_symbolic_model(abstract_system)
println(SY.is_deterministic(deterministic_abstract_system))

using JuMP
new_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("abstract_system"), deterministic_abstract_system)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(new_optimizer, MOI.RawOptimizerAttribute("sparse_input"), true)
MOI.optimize!(new_optimizer)

concrete_controller = MOI.get(new_optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_problem = MOI.get(new_optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
x0 = SVector(0.4, 0.4, 0.0)
new_control_trajectory = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

fig = plot(; aspect_ratio = :equal)
plot!(concrete_system.X; color = :yellow, opacity = 0.5, label = "")
plot!(deterministic_abstract_system.Xdom; color = :blue, opacity = 0.5, label = "", efficient=false)
plot!(concrete_problem.initial_set; color = :green, label = "", opacity = 0.2)
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, label = "", opacity = 0.2)
plot!(
        SY.get_domain_from_states(deterministic_abstract_system, abstract_problem.initial_set);
        color = :green,
        label = "",
    )
plot!(
        SY.get_domain_from_states(deterministic_abstract_system, abstract_problem.target_set);
        color = :red,
        label = "",
    )
plot!(control_trajectory; ms = 0.5)
plot!(new_control_trajectory; ms = 0.5, color = :red)
display(fig)

end


