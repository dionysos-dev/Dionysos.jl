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
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/non_linear.jl")

U = UT.HyperRectangle(SVector(-10.0, -10.0), SVector(10.0, 10.0))
xnew = SVector{2, Float64}([1.0; 1.0])
ρ = 0.00005
Wbound = 0.0
λ = 0.01

concrete_problem = NonLinear.problem(;
    X = IntervalBox(-20.0 .. 20.0, 2),
    obstacles = [
        UT.Ellipsoid(Matrix{Float64}(I(2)) * 1 / 50, [0.0; 0.0]),
        UT.Ellipsoid([0.2 0.2; 0.2 2.0] * 0.4, [15.0; -7.0]),
        UT.Ellipsoid([2.0 0.2; 0.2 0.5] * 0.2, [20.0; 0.0]),
    ],
    U = U,
    E0 = UT.Ellipsoid(Matrix{Float64}(I(2)) * 10.0, [-10.0; -10.0]),
    Ef = UT.Ellipsoid(Matrix{Float64}(I(2)) * 1.0, [10.0; 10.0]),
    state_cost = UT.ZeroFunction(),
    transition_cost = UT.QuadraticStateControlFunction(
        Matrix{Float64}(I(2)),
        Matrix{Float64}(I(2)),
        zeros(2, 2),
        zeros(2),
        zeros(2),
        1.0,
    ),
    W = UT.HyperRectangle(SVector(-Wbound, -Wbound), SVector(Wbound, Wbound)),
    noise = false,
    μ = ρ,
)

concrete_system = concrete_problem.system

# Optimizer's parameters
const FALLBACK_URL = "mosek://solve.mosek.com:30080"
sdp_opt = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
MOI.set(sdp_opt, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)
maxδx = 100
maxδu = 10 * 2
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
start_time = time()
MOI.optimize!(optimizer)
elapsed_time = time() - start_time
println("Elapsed time: ", elapsed_time, " seconds")

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
nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time;
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
cost_true = sum(cost_control_trajectory.costs.seq);
println("Goal set reached")
println("Guaranteed cost:\t $(cost_bound)")
println("True cost:\t\t $(cost_true)")

# ## Display the results
# # Display the specifications and domains
fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
)

#Display the concrete domain
for obs in concrete_system.obstacles
    plot!(obs; color = :black)
end

#Display the abstract domain
plot!(concrete_problem.target_set; color = :red)
plot!(abstract_system; arrowsB = true, cost = true)

#Display the concrete specifications
plot!(concrete_problem.initial_set; color = :green)

xlabel!("\$x_1\$")
ylabel!("\$x_2\$")
xlims!(-25, 25)

display(fig)
