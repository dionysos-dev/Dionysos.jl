using StaticArrays, LinearAlgebra, Plots
using JuMP, Clarabel
import CDDLib

import Random
Random.seed!(0)

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

lib = CDDLib.Library() # polyhedron lib

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "pwa_sys.jl"))

Usz = 70 # upper limit on |u|, `Usz = 50` in [1]
Wsz = 3 # `Wsz = 5` in [1]
dt = 0.01; # discretization step

concrete_problem =
    PWAsys.problem(; lib = lib, dt = dt, Usz = Usz, Wsz = Wsz, simple = false)
concrete_system = concrete_problem.system

n_step = 3
origin = SVector(0.0, 0.0)
h = SVector(1.0 / n_step, 1.0 / n_step)
nx = size(concrete_system.resetmaps[1].A, 1)

P = (1 / nx) * diagm((h ./ 2) .^ (-2))
state_grid = MP.GridEllipsoidalRectangular(origin, h, P)

R = h ./ 2

Pm = P

opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

optimizer = MOI.instantiate(AB.UniformEllipsoidAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
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
);

MOI.optimize!(optimizer)

abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
concrete_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_value_function"))
abstract_value_function =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_value_function"));

get_mode(x) = findfirst(m -> (x ∈ m.X), concrete_system.resetmaps)

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

nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time;
function reached(x)
    states = SY.get_abstract_states(abstract_system, x)
    if !isempty(states ∩ abstract_problem.target_set)
        return true
    else
        return false
    end
end

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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
