using StaticArrays, LinearAlgebra, Random, IntervalArithmetic
using MathematicalSystems, HybridSystems
using JuMP, Mosek, MosekTools
using Plots, Colors
using Test
Random.seed!(0)

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const CO = DI.Control
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include("../../../problems/non_linear.jl")

concrete_problem = NonLinear.problem()
concrete_system = concrete_problem.system

const FALLBACK_URL = "mosek://solve.mosek.com:30080"
sdp_opt = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
MOI.set(sdp_opt, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)

maxδx = 100 # 100
maxδu = 10 * 2 # Usz * 2
λ = 0.01 # 0.01
k1 = 1
k2 = 1
RRTstar = false
continues = false
maxIter = 100 # 100

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

using Suppressor
@suppress begin # this is a workaround to supress the undesired output of SDPA
    MOI.optimize!(optimizer)
end

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_lyap_fun"))
concrete_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_lyap_fun"));

cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost, x, u)
reached(x) = x ∈ concrete_problem.target_set
nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time; # max num of steps

x0 = concrete_problem.initial_set.c
x_traj, u_traj, cost_traj = CO.get_closed_loop_trajectory(
    concrete_system.f_eval,
    concrete_controller,
    cost_eval,
    x0,
    nstep;
    stopping = reached,
    noise = true,
)
cost_bound = concrete_lyap_fun(x0)
cost_true = sum(cost_traj);
println("Goal set reached")
println("Guaranteed cost:\t $(cost_bound)")
println("True cost:\t\t $(cost_true)")

fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
);
xlabel!("\$x_1\$");
ylabel!("\$x_2\$");
title!("Specifictions and domains");

#Display the concrete domain
plot!(concrete_system.X; color = :yellow, opacity = 0.5);
for obs in concrete_system.obstacles
    plot!(obs; color = :black)
end

#Display the abstract domain
plot!(abstract_system; arrowsB = false, cost = false);

#Display the concrete specifications
plot!(concrete_problem.initial_set; color = :green);
plot!(concrete_problem.target_set; color = :red)

fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
);
title!("Abstractions");
plot!(abstract_system; arrowsB = true, cost = false)

fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
);
xlabel!("\$x_1\$");
ylabel!("\$x_2\$");
title!("Trajectory and Lyapunov-like Fun.");

for obs in concrete_system.obstacles
    plot!(obs; color = :black)
end
plot!(abstract_system; arrowsB = false, cost = true);
plot!(UT.DrawTrajectory(x_traj); color = :black)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
