using StaticArrays, LinearAlgebra, Random, IntervalArithmetic
using MathematicalSystems, HybridSystems
using JuMP, Clarabel
using Plots, Colors
using Test
Random.seed!(0)

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include("../../../problems/non_linear.jl")

concrete_problem = NonLinear.problem()
concrete_system = concrete_problem.system

sdp_opt = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

maxδx = 100
maxδu = 10 * 2
λ = 0.01
k1 = 1
k2 = 1
RRTstar = false
continues = false
maxIter = 100

optimizer = MOI.instantiate(AB.LazyEllipsoidsAbstraction.Optimizer)
AB.LazyEllipsoidsAbstraction.set_optimizer!(
    optimizer,
    concrete_problem,
    sdp_opt,
    maxδx,
    maxδu,
    λ,
    k1,
    k2,
    RRTstar,
    continues,
    maxIter,
)

# Build the state feedback abstraction and solve the optimal control problem using RRT algorithm.
MOI.optimize!(optimizer)

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_lyap_fun"))
concrete_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_lyap_fun"));

# ## Simulation
# We define the cost and stopping criteria for a simulation
cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost, x, u)
reached(x) = x ∈ concrete_problem.target_set
nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time; # max num of steps
# We simulate the closed loop trajectory
x0 = concrete_problem.initial_set.c
cost_control_trajectory = ST.get_closed_loop_trajectory(
    concrete_system.f_eval,
    concrete_controller,
    cost_eval,
    x0,
    nstep;
    stopping = reached,
    noise = true,
)
cost_bound = concrete_lyap_fun(x0)
cost_true = ST.get_cost(cost_control_trajectory)
println("Goal set reached")
println("Guaranteed cost:\t $(cost_bound)")
println("True cost:\t\t $(cost_true)")

# ## Display the results
# # Display the specifications and domains
fig1 = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
)
xlabel!("\$x_1\$")
ylabel!("\$x_2\$")
title!("Specifictions and domains")

#Display the concrete domain
plot!(fig1, concrete_system.X; color = :yellow, opacity = 0.5)
for obs in concrete_system.obstacles
    plot!(fig1, obs; color = :black)
end

#Display the abstract domain
plot!(fig1, abstract_system; arrowsB = false, cost = false)

#Display the concrete specifications
plot!(fig1, concrete_problem.initial_set; color = :green)
plot!(fig1, concrete_problem.target_set; color = :red)

# # Display the abstraction
fig2 = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
)
title!("Abstractions")
plot!(fig2, abstract_system; arrowsB = true, cost = false)

# # Display the Lyapunov function and the trajectory
fig3 = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
)
xlabel!("\$x_1\$")
ylabel!("\$x_2\$")
title!("Trajectory and Lyapunov-like Fun.")

for obs in concrete_system.obstacles
    plot!(fig3, obs; color = :black)
end
plot!(fig3, abstract_system; arrowsB = false, cost = true)
plot!(fig3, cost_control_trajectory; color = :black)
