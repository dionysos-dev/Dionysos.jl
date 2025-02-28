module TestMain
using Test

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

# Example of ellipsoidal based abstraction

if !isdefined(@__MODULE__, :Usz)
    Usz = 70 # upper limit on |u|
    Wsz = 3
    dt = 0.01 # discretization step
    n_step = 3 # discretization of one unit of space
    simple = true
    no_plot = true
end
lib = CDDLib.Library() # polyhedron lib
include("../../problems/pwa_sys.jl")

# Problem parameters

# Usz = 50 # upper limit on |u|
# Wsz = 5
# n_step = 5 # discretization of one unit of space
# simple = false
# no_plot = false

# Usz = 50 # upper limit on |u|
# Wsz = 5
# dt = 0.01; # discretization step
# n_step = 5 # discretization of one unit of space
# simple = false
# no_plot = false

opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

concrete_problem =
    PWAsys.problem(; lib = lib, dt = dt, Usz = Usz, Wsz = Wsz, simple = simple)
concrete_system = concrete_problem.system

# Abstraction parameters
X_origin = SVector(0.0, 0.0)
X_step = SVector(1.0 / n_step, 1.0 / n_step)
nx = size(concrete_system.resetmaps[1].A, 1)
P = (1 / nx) * diagm((X_step ./ 2) .^ (-2))
state_grid = DO.GridEllipsoidalRectangular(X_origin, X_step, P)

optimizer = MOI.instantiate(AB.EllipsoidsAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("sdp_solver"), opt_sdp)

# Build the state-feedback abstraction and solve the optimal control problem by through Dijkstra's algorithm [2, p.86].
MOI.optimize!(optimizer)

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_lyap_fun"))
concrete_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_lyap_fun"))
transitionCont = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionCont"))
transitionCost = MOI.get(optimizer, MOI.RawOptimizerAttribute("transitionCost"))

# return pwa mode for a given x
get_mode(x) = findfirst(m -> (x ∈ m.X), concrete_system.resetmaps)
function f_eval1(x, u)
    currState = SY.get_states_by_xpos(
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

### Simulation
# We define the stopping criteria for a simulation
nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time #max num of steps
function reached(x)
    currState = SY.get_states_by_xpos(
        abstract_system,
        DO.crop_to_domain(abstract_system.Xdom, DO.get_all_pos_by_coord(state_grid, x)),
    )
    if !isempty(currState ∩ abstract_problem.target_set)
        return true
    else
        return false
    end
end

# We simulate the closed loop trajectory
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
cost_true = sum(cost_control_trajectory.costs.seq)
println("Goal set reached")
println("Guaranteed cost:\t $(cost_bound)")
println("True cost:\t\t $(cost_true)")

@static if get(ENV, "CI", "false") == "false" &&
           (isdefined(@__MODULE__, :no_plot) && no_plot == false)
    ## Display the specifications and domains
    fig = plot(;
        aspect_ratio = :equal,
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 16,
        titlefontsize = 14,
    )
    rectX = concrete_system.ext[:X]
    xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2)
    ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2)

    # We display the concrete domain
    plot!(rectX; color = :yellow, opacity = 0.5)

    # We display the abstract domain
    plot!(abstract_system.Xdom; color = :blue, opacity = 0.5)

    # We display the abstract specifications
    plot!(
        SY.get_domain_from_symbols(abstract_system, abstract_problem.initial_set);
        color = :green,
        opacity = 0.5,
    )
    plot!(
        SY.get_domain_from_symbols(abstract_system, abstract_problem.target_set);
        color = :red,
        opacity = 0.5,
    )

    # We display the concrete specifications
    plot!(UT.DrawPoint(concrete_problem.initial_set); color = :green, opacity = 1.0)
    plot!(UT.DrawPoint(concrete_problem.target_set); color = :red, opacity = 1.0)

    xlabel!("\$x_1\$")
    ylabel!("\$x_2\$")
    title!("Specifictions and domains")
    display(fig)

    ## Display the abstraction
    fig = plot(;
        aspect_ratio = :equal,
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 16,
        titlefontsize = 14,
    )
    xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2)
    ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2)
    plot!(abstract_system; arrowsB = true, cost = false)
    title!("Abstractions")
    display(fig)

    ## Display the Lyapunov function and the trajectory
    fig = plot(;
        aspect_ratio = :equal,
        xtickfontsize = 10,
        ytickfontsize = 10,
        guidefontsize = 16,
        titlefontsize = 14,
    )
    xlims!(rectX.A.lb[1] - 0.2, rectX.A.ub[1] + 0.2)
    ylims!(rectX.A.lb[2] - 0.2, rectX.A.ub[2] + 0.2)
    plot!(abstract_system; arrowsB = false, cost = true, lyap_fun = optimizer.lyap)
    plot!(cost_control_trajectory; color = :black)
    xlabel!("\$x_1\$")
    ylabel!("\$x_2\$")
    title!("Trajectory and Lyapunov-like Fun.")
    display(fig)
end
@testset "state_trans" begin
    @test cost_bound ≈ 0.6250139513432214 rtol = 1e-3
    @test cost_true ≈ 0.36844089806471475 rtol = 1e-1
    @test cost_true <= cost_bound
end
end
