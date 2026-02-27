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
include("../problems/pwa_sys.jl")

Usz = 70
Wsz = 3
dt = 0.01

concrete_problem =
    PWAsys.problem(; lib = lib, dt = dt, Usz = Usz, Wsz = Wsz, simple = true)

concrete_system = concrete_problem.system

# ---------------------------------------------------------
# Build EmptyProblem (abstraction-only)
# ---------------------------------------------------------

empty_problem = PR.EmptyProblem(
    concrete_system,
    concrete_system.ext[:X]   # region = state constraint set
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
opt_sdp = optimizer_with_attributes(
    Clarabel.Optimizer,
    MOI.Silent() => true
)

# ---------------------------------------------------------
# Instantiate abstraction optimizer
# ---------------------------------------------------------

optimizer = MOI.instantiate(AB.UniformEllipsoidAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("empty_problem"), empty_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("P"), P)
MOI.set(optimizer, MOI.RawOptimizerAttribute("Pm"), Pm)
MOI.set(optimizer, MOI.RawOptimizerAttribute("R"), R)
MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), opt_sdp)

# ---------------------------------------------------------
# Build abstraction
# ---------------------------------------------------------

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstraction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")

println("Abstraction built.")
println("Number of states: ", SY.get_n_state(abstract_system))
println("Number of transitions: ", length(optimizer.abstraction_solver.transitionCost))

XMapping = SY.get_state_mapping(abstract_system)
Xset = SY.get_state_domain(abstract_system)

fig = plot(; aspect_ratio = :equal);
# plot!(XMapping; efficient=true, color=:grey)
# plot!((Xset, XMapping); efficient=false, color=:grey)
plot!(abstract_system; efficient=false, with_arrows=true)
display(fig)

# ---------------------------------------------------------
# Solve control problem
# ---------------------------------------------------------

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

MOI.optimize!(optimizer)
abstract_problem_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")


# # ## Define the mapping function
# #Return pwa mode for a given x
# get_mode(x) = findfirst(m -> (x ∈ m.X), concrete_system.resetmaps)
# # To simplify : "We assume that inside cells intersecting the boundary of partitions of X the selected piecewise-affine mode is the same all over its interior and given by the mode
# # defined at its center."
# function f_eval1(x, u)
#     states = SY.get_states_by_xpos(
#         abstract_system,
#         DO.crop_to_domain(abstract_system.Xdom, DO.get_all_pos_by_coord(state_grid, x)),
#     )
#     from = nothing
#     for s in states
#         if s in abstract_controller.X
#             from = s
#             break
#         end
#     end
#     c = DO.get_coord_by_pos(state_grid, SY.get_xpos_by_state(abstract_system, from))
#     m = get_mode(c)
#     W = concrete_system.ext[:W]
#     w = (2 * (rand(2) .^ (1 / 4)) .- 1) .* W[:, 1]
#     return concrete_system.resetmaps[m].A * x +
#            concrete_system.resetmaps[m].B * u +
#            concrete_system.resetmaps[m].c +
#            w
# end

# cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost[1][1], x, u)

# # ### Simulation

# # We define the stopping criteria for a simulation
# nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time; #max num of steps
# function reached(x)
#     states = SY.get_states_by_xpos(
#         abstract_system,
#         DO.crop_to_domain(abstract_system.Xdom, DO.get_all_pos_by_coord(state_grid, x)),
#     )
#     if !isempty(states ∩ abstract_problem.target_set)
#         return true
#     else
#         return false
#     end
# end

# # We simulate the closed loop trajectory
# x0 = concrete_problem.initial_set
# x_traj, u_traj = ST.get_closed_loop_trajectory(
#     concrete_system,
#     concrete_controller,
#     x0,
#     nstep;
#     stopping = reached,
#     f_map_override = f_eval1,
# )
# c_traj, cost_true = ST.get_cost_trajectory(x_traj, u_traj, cost_eval)
# cost_bound = concrete_lyap_fun(x0)
# println("Goal set reached")
# println("Guaranteed cost:\t $(cost_bound)")
# println("True cost:\t\t $(cost_true)")

# # ### Visualize the results. 
# rectX = concrete_system.ext[:X];

# # ## Display the specifications and domains
# fig = plot(;
#     aspect_ratio = :equal,
#     xtickfontsize = 10,
#     ytickfontsize = 10,
#     guidefontsize = 16,
#     titlefontsize = 14,
# );
# xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2);
# ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2);
# xlabel!("\$x_1\$");
# ylabel!("\$x_2\$");
# title!("Specifictions and domains");
# #We display the concrete domain
# plot!(rectX; color = :grey, opacity = 1.0, label = "");
# #We display the abstract domain
# plot!(abstract_system.Xdom; color = :blue, efficient = false, opacity = 0.5);
# #We display the abstract specifications
# plot!(
#     SY.get_domain_from_states(abstract_system, abstract_problem.initial_set);
#     color = :green,
#     efficient = false,
#     opacity = 0.5,
# );
# plot!(
#     SY.get_domain_from_states(abstract_system, abstract_problem.target_set);
#     color = :red,
#     efficient = false,
#     opacity = 0.5,
# );
# #We display the concrete specifications
# plot!(UT.DrawPoint(concrete_problem.initial_set); color = :green, opacity = 1.0);
# plot!(UT.DrawPoint(concrete_problem.target_set); color = :red, opacity = 1.0)

# # ## Display the abstraction
# fig = plot(;
#     aspect_ratio = :equal,
#     xtickfontsize = 10,
#     ytickfontsize = 10,
#     guidefontsize = 16,
#     titlefontsize = 14,
# );
# xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2);
# ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2);
# title!("Abstractions");
# plot!(abstract_system; arrowsB = true, efficient = false, cost = false)

# # ## Display the Lyapunov function and the trajectory
# fig = plot(;
#     aspect_ratio = :equal,
#     xtickfontsize = 10,
#     ytickfontsize = 10,
#     guidefontsize = 16,
#     titlefontsize = 14,
# );
# xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2);
# ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2);
# xlabel!("\$x_1\$");
# ylabel!("\$x_2\$");
# title!("Trajectory and Lyapunov-like Fun.");
# plot!(abstract_system; arrowsB = false, value_function = optimizer.abstract_lyap_fun);
# plot!(x_traj; color = :black)

# # ## References
# #
# # 1. L. N. Egidio, T. Alves Lima, R. M. Jungers, "State-feedback Abstractions for Optimal Control of Piecewise-affine Systems", IEEE 61st Conference on Decision and Control (CDC), 2022, accepted.
# # 1. D. Bertsekas, "Dynamic programming and optimal control". Volume I, Athena scientific, 2012.
