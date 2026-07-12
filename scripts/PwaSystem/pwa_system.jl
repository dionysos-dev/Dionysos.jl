using StaticArrays, LinearAlgebra
using JuMP, Clarabel
import MathOptInterface as MOI
import CDDLib
import Random
Random.seed!(0)

using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

# ---------------------------------------------------------
# Load PWA system
# ---------------------------------------------------------

lib = CDDLib.Library()
include("../../problems/PwaSystem/pwa_system.jl")

Usz = 70
Wsz = 3
dt = 0.01

concrete_problem =
    PwaSystem.problem(; lib = lib, dt = dt, Usz = Usz, Wsz = Wsz, simple = false)

concrete_system = concrete_problem.system

# ---------------------------------------------------------
# Build AlternatingSimulationProblem (abstraction-only)
# ---------------------------------------------------------

alternating_simulation_problem = PR.AlternatingSimulationProblem(
    concrete_system,
    concrete_system.ext[:X],   # region = state constraint set
)

# ---------------------------------------------------------
# Ellipsoidal Grid
# ---------------------------------------------------------

n_step = 3

origin = SVector(0.0, 0.0)
h = SVector(1.0 / n_step, 1.0 / n_step)

nx = size(concrete_system.resetmaps[1].A, 1)

# Ellipsoid matrix P
P = (1 / nx) * diagm((h ./ 2) .^ (-2))

state_grid = MP.GridEllipsoidalRectangular(origin, h, P)

# Overapproximation radius (same size as state)
R = h ./ 2

# Optional: same shape for source and target ellipsoids
Pm = P

# SDP solver
opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

# ---------------------------------------------------------
# Instantiate abstraction optimizer
# ---------------------------------------------------------

optimizer = MOI.instantiate(AB.UniformEllipsoidAbstraction.Optimizer)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("alternating_simulation_problem"),
    alternating_simulation_problem,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("incl_mode"), MP.INNER)
MOI.set(optimizer, MOI.RawOptimizerAttribute("P"), P)
MOI.set(optimizer, MOI.RawOptimizerAttribute("Pm"), Pm)
MOI.set(optimizer, MOI.RawOptimizerAttribute("R"), R)
MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), opt_sdp)
nx, nu = 2, 2
naug = nx + nu + 1
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("Q_aug"),
    Matrix{Float64}(LinearAlgebra.I, naug, naug)*(dt^2),
)

# ---------------------------------------------------------
# Build abstraction
# ---------------------------------------------------------

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")

println("Abstraction built.")
println("Number of states: ", SY.get_n_state(abstract_system))
println("Number of transitions: ", length(optimizer.abstraction_solver.transitionCost))

XMapping = SY.get_state_mapping(abstract_system)
Xset = SY.get_state_set(abstract_system)

fig = plot(; aspect_ratio = :equal);
# plot!(XMapping; efficient=true, color=:grey)
# plot!((Xset, XMapping); efficient=false, color=:grey)
plot!(abstract_system; efficient = false, with_arrows = true)
display(fig)

# ---------------------------------------------------------
# Solve control problem
# ---------------------------------------------------------

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

MOI.optimize!(optimizer)
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
concrete_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_value_function"))
abstract_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"))
println("Time to solve the abstract problem: $(abstract_problem_time)")

# # ### Simulation

# Return pwa mode for a given x
get_mode(x) = findfirst(m -> (x ∈ m.X), concrete_system.resetmaps)
# To simplify : "We assume that inside cells intersecting the boundary of partitions of X the selected piecewise-affine mode is the same all over its interior and given by the mode
# defined at its center."
function f_eval1(x, u)
    states = SY.get_abstract_states(abstract_system, x)
    min_state = argmin(s -> abstract_value_function(s), states)
    c = SY.get_concrete_state(abstract_system, min_state)
    m = get_mode(c)
    W = concrete_system.ext[:W]
    w = (2 * (rand(2) .^ (1 / 4)) .- 1) .* W[:, 1]
    return concrete_system.resetmaps[m].A * x +
           concrete_system.resetmaps[m].B * u +
           concrete_system.resetmaps[m].c +
           w
end
cost_eval(x, u) = concrete_problem.transition_cost[1][1](x, u)

# We define the stopping criteria for a simulation
nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time;
function reached(x)
    states = SY.get_abstract_states(abstract_system, x)
    if !isempty(states ∩ abstract_problem.target_set)
        return true
    else
        return false
    end
end

# We simulate the closed loop trajectory
x0 = concrete_problem.initial_set
x_traj, u_traj = ST.get_closed_loop_trajectory(
    concrete_system,
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
    f_map_override = f_eval1,
)
c_traj, cost_true = ST.get_cost_trajectory(x_traj, u_traj, cost_eval)
cost_bound = concrete_value_function(x0)
println("Goal set reached")
println("Guaranteed cost:\t $(cost_bound)")
println("True cost:\t\t $(cost_true)")

# ### Visualize the results. 

Xmap = SY.get_state_mapping(abstract_system)
fig = plot(; aspect_ratio = :equal)
X = concrete_system.ext[:X]
plot!(X; color = :grey, opacity = 1.0, label = "")
plot!(abstract_system; value_function = abstract_value_function)
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.initial_set), Xmap);
    color = :green,
    efficient = false,
    opacity = 0.6,
)
plot!(
    (SY.get_state_set_from_states(abstract_system, abstract_problem.target_set), Xmap);
    color = :red,
    efficient = false,
    opacity = 0.6,
)
plot!(UT.DrawPoint(concrete_problem.initial_set); color = :green, opacity = 1.0);
plot!(UT.DrawPoint(concrete_problem.target_set); color = :red, opacity = 1.0)
plot!(x_traj; ms = 2.0, arrows = false, color = :blue)
display(fig)
