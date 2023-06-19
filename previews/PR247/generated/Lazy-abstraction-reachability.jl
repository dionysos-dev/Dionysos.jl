using StaticArrays, JuMP, Plots

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

include("../../../problems/simple_problem.jl")

# specific functions
function post_image(abstract_system, concrete_system, xpos, u)
    Xdom = abstract_system.Xdom
    x = DO.get_coord_by_pos(Xdom.grid, xpos)
    Fx = concrete_system.f_eval(x, u)
    r = Xdom.grid.h / 2.0 + concrete_system.measnoise
    Fr = r

    rectI = DO.get_pos_lims_outer(Xdom.grid, UT.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(DO._ranges(rectI)...)
    over_approx = []
    allin = true
    for ypos in ypos_iter
        ypos = DO.set_in_period_pos(Xdom, ypos)
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = SY.get_state_by_xpos(abstract_system, ypos)[1]
        push!(over_approx, target)
    end
    return allin ? over_approx : []
end

function pre_image(abstract_system, concrete_system, xpos, u)
    grid = abstract_system.Xdom.grid
    x = DO.get_coord_by_pos(grid, xpos)
    potential = Int[]
    x_prev = concrete_system.f_backward(x, u)
    xpos_cell = DO.get_pos_by_coord(grid, x_prev)
    n = 2
    for i in (-n):n
        for j in (-n):n
            x_n = (xpos_cell[1] + i, xpos_cell[2] + j)
            x_n = DO.set_in_period_pos(abstract_system.Xdom, x_n)
            if x_n in abstract_system.Xdom
                cell = SY.get_state_by_xpos(abstract_system, x_n)[1]
                if !(cell in potential)
                    push!(potential, cell)
                end
            end
        end
    end
    return potential
end

function compute_reachable_set(rect::UT.HyperRectangle, concrete_system, Udom)
    r = (rect.ub - rect.lb) / 2.0 + concrete_system.measnoise
    Fr = r
    x = UT.get_center(rect)
    n = UT.get_dims(rect)
    lb = fill(Inf, n)
    ub = fill(-Inf, n)
    for upos in DO.enum_pos(Udom)
        u = DO.get_coord_by_pos(Udom.grid, upos)
        Fx = concrete_system.f_eval(x, u)
        lb = min.(lb, Fx .- Fr)
        ub = max.(ub, Fx .+ Fr)
    end
    lb = SVector{n}(lb)
    ub = SVector{n}(ub)
    return UT.HyperRectangle(lb, ub)
end

concrete_problem = SimpleProblem.problem()
concrete_system = concrete_problem.system

hx = [0.5, 0.5]
u0 = SVector(0.0, 0.0)
hu = SVector(0.5, 0.5)
Ugrid = DO.GridFree(u0, hu)
maxIter = 100

optimizer = MOI.instantiate(AB.LazyAbstraction.Optimizer)

AB.LazyAbstraction.set_optimizer!(
    optimizer,
    concrete_problem,
    maxIter,
    pre_image,
    post_image,
    compute_reachable_set,
    hx,
    Ugrid,
)

using Suppressor
@suppress begin # this is a workaround to supress the undesired output of SDPA
    MOI.optimize!(optimizer)
end

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost, x, u)
reached(x) = x ∈ concrete_problem.target_set
nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time; # max num of steps

x0 = UT.get_center(concrete_problem.initial_set)
x_traj, u_traj, cost_traj = CO.get_closed_loop_trajectory(
    concrete_system.f_eval,
    concrete_controller,
    cost_eval,
    x0,
    nstep;
    stopping = reached,
    noise = false,
)

cost_true = sum(cost_traj);
println("Goal set reached")
println("True cost:\t\t $(cost_true)")

fig = plot(; aspect_ratio = :equal);

plot!(concrete_system.X; color = :yellow, opacity = 0.5);

plot!(abstract_system.Xdom; color = :blue, opacity = 0.5);

plot!(concrete_problem.initial_set; color = :green, opacity = 0.8);
plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.8);

plot!(UT.DrawTrajectory(x_traj); ms = 0.5)

fig = plot(; aspect_ratio = :equal)
x0 = SVector(5.5, 5.5)
AB.LazyAbstraction.plot_result!(optimizer.lazy_search_problem; x0 = x0)
plot!(; show = true, legend = false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

