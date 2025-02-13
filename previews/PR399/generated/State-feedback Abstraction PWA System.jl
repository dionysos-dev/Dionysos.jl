using StaticArrays, LinearAlgebra, Polyhedra, Random
using MathematicalSystems, HybridSystems
using JuMP, Clarabel
using SemialgebraicSets, CDDLib
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

lib = CDDLib.Library() # polyhedron lib

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "pwa_sys.jl"))

Usz = 70 # upper limit on |u|, `Usz = 50` in [1]
Wsz = 3 # `Wsz = 5` in [1]
dt = 0.01; # discretization step

concrete_problem =
    PWAsys.problem(; lib = lib, dt = dt, Usz = Usz, Wsz = Wsz, simple = false)
concrete_system = concrete_problem.system

n_step = 3
X_origin = SVector(0.0, 0.0);
X_step = SVector(1.0 / n_step, 1.0 / n_step)
nx = size(concrete_system.resetmaps[1].A, 1)
P = (1 / nx) * diagm((X_step ./ 2) .^ (-2))
state_grid = DO.GridEllipsoidalRectangular(X_origin, X_step, P);
opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

optimizer = MOI.instantiate(AB.EllipsoidsAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), opt_sdp)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_lyap_fun"))
concrete_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_lyap_fun"))
transitionCont = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionCont"))
transitionCost = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionCost"));

#Return pwa mode for a given x
get_mode(x) = findfirst(m -> (x ∈ m.X), concrete_system.resetmaps)

function f_eval1(x, u)
    currState = SY.get_all_states_by_xpos(
        abstract_system,
        DO.crop_to_domain(abstract_system.Xdom, DO.get_all_pos_by_coord(state_grid, x)),
    )
    next_action = nothing
    for action in abstract_controller.data
        if (action[1] ∩ currState) ≠ []
            next_action = action
        end
    end
    c = DO.get_coord_by_pos(
        state_grid,
        SY.get_xpos_by_state(abstract_system, next_action[1]),
    )
    m = get_mode(c)
    W = concrete_system.ext[:W]
    w = (2 * (rand(2) .^ (1 / 4)) .- 1) .* W[:, 1]
    return concrete_system.resetmaps[m].A * x +
           concrete_system.resetmaps[m].B * u +
           concrete_system.resetmaps[m].c +
           w
end

cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost[1][1], x, u)

nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time; #max num of steps
function reached(x)
    currState = SY.get_all_states_by_xpos(
        abstract_system,
        DO.crop_to_domain(abstract_system.Xdom, DO.get_all_pos_by_coord(state_grid, x)),
    )
    if !isempty(currState ∩ abstract_problem.target_set)
        return true
    else
        return false
    end
end

x0 = concrete_problem.initial_set
cost_control_trajectory = ST.get_closed_loop_trajectory(
    f_eval1,
    concrete_controller,
    cost_eval,
    x0,
    nstep;
    stopping = reached,
)
cost_bound = concrete_lyap_fun(x0)
cost_true = sum(cost_control_trajectory.costs.seq);
println("Goal set reached")
println("Guaranteed cost:\t $(cost_bound)")
println("True cost:\t\t $(cost_true)")

rectX = concrete_system.ext[:X];

fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
);
xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2);
ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2);
xlabel!("\$x_1\$");
ylabel!("\$x_2\$");
title!("Specifictions and domains");
#We display the concrete domain
plot!(rectX; color = :yellow, opacity = 0.5);
#We display the abstract domain
plot!(abstract_system.Xdom; color = :blue, opacity = 0.5);
#We display the abstract specifications
plot!(
    SY.get_domain_from_states(abstract_system, abstract_problem.initial_set);
    color = :green,
    opacity = 0.5,
);
plot!(
    SY.get_domain_from_states(abstract_system, abstract_problem.target_set);
    color = :red,
    opacity = 0.5,
);
#We display the concrete specifications
plot!(UT.DrawPoint(concrete_problem.initial_set); color = :green, opacity = 1.0);
plot!(UT.DrawPoint(concrete_problem.target_set); color = :red, opacity = 1.0)

fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
);
xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2);
ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2);
title!("Abstractions");
plot!(abstract_system; arrowsB = true, cost = false)

fig = plot(;
    aspect_ratio = :equal,
    xtickfontsize = 10,
    ytickfontsize = 10,
    guidefontsize = 16,
    titlefontsize = 14,
);
xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2);
ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2);
xlabel!("\$x_1\$");
ylabel!("\$x_2\$");
title!("Trajectory and Lyapunov-like Fun.");
plot!(abstract_system; arrowsB = false, cost = true, lyap_fun = optimizer.lyap);
plot!(cost_control_trajectory; color = :black)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
