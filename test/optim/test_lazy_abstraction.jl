module TestMain
using StaticArrays, JuMP, Plots, MathematicalSystems
MS = MathematicalSystems
using Test
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

rectX = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
obsX = UT.LazyUnionSetArray([UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))])
_X_ = UT.LazySetMinus(rectX, obsX)

rectU = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
obsU = UT.LazyUnionSetArray([UT.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))])
_U_ = UT.LazySetMinus(rectU, obsU)

dt = 0.8
fd = (x, u) -> x + dt*u
concrete_system = MS.ConstrainedBlackBoxControlDiscreteSystem(fd, 2, 2, _X_, _U_)

_I_ = UT.HyperRectangle(SVector(5.0, 5.0), SVector(6.0, 6.0))
_T_ = UT.HyperRectangle(SVector(25.0, 25.0), SVector(28.0, 28.0))
concrete_problem = PR.OptimalControlProblem(
    concrete_system,
    _I_,
    _T_,
    UT.ZeroFunction(),
    UT.ConstantControlFunction(0.5),
    PR.Infinity(),
)

function post_image(abstract_system, concrete_system, xpos, u)
    f = MS.mapping(concrete_system)   # (x,u) -> xnext

    Xdom = abstract_system.Xdom
    x = DO.get_coord_by_pos(Xdom, xpos)
    Fx = f(x, u)
    r = DO.get_grid(Xdom).h / 2.0
    Fr = r

    rectI = DO.get_pos_lims_outer(DO.get_grid(Xdom), UT.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(DO._ranges(rectI)...)
    over_approx = []
    allin = true
    for ypos in ypos_iter
        if !(ypos in abstract_system)
            allin = false
            break
        end
        target = SY.get_state_by_xpos(abstract_system, ypos)[1]
        push!(over_approx, target)
    end
    return allin ? over_approx : []
end

## specific functions
backward_map = (x, u) -> x - dt*u
function pre_image(abstract_system, xpos, u)
    Xdom = abstract_system.Xdom
    x = DO.get_coord_by_pos(Xdom, xpos)

    potential = Int[]
    x_prev = backward_map(x, u)
    xpos_cell = DO.get_pos_by_coord(Xdom, x_prev)

    n = 2
    for i in (-n):n, j in (-n):n
        x_n = (xpos_cell[1] + i, xpos_cell[2] + j)
        if x_n in abstract_system
            cell = SY.get_state_by_xpos(abstract_system, x_n)[1]
            if !(cell in potential)
                push!(potential, cell)
            end
        end
    end
    return potential
end

function compute_reachable_set(rect::UT.HyperRectangle, concrete_system, Udom)
    f = MS.mapping(concrete_system)   # (x,u) -> xnext

    r = (rect.ub - rect.lb) / 2.0
    Fr = r
    x = UT.get_center(rect)
    n = UT.get_dims(rect)

    lb = fill(Inf, n)
    ub = fill(-Inf, n)

    for u in DO.enum_elems(Udom)
        Fx = f(x, u)
        lb = min.(lb, Fx .- Fr)
        ub = max.(ub, Fx .+ Fr)
    end
    return UT.HyperRectangle(SVector{n}(lb), SVector{n}(ub))
end
minimum_transition_cost(symmodel, contsys, source, target) = 1.0

concrete_system = concrete_problem.system

hx = SVector(0.5, 0.5)
u0 = SVector(0.0, 0.0)
hu = SVector(0.5, 0.5)
Ugrid = DO.GridFree(u0, hu)
hx_heuristic = SVector(1.5, 1.5)
maxIter = 100
periodic_dims = SVector{0, Int}()
periodic_periods = SVector{0, Float64}()
periodic_start = SVector{0, Float64}()

# setup for the tests
if !isdefined(@__MODULE__, :hx)
    hx = SVector(0.5, 0.5)
    u0 = SVector(0.0, 0.0)
    hu = SVector(0.5, 0.5)
    Ugrid = DO.GridFree(u0, hu)
    hx_heuristic = SVector(1.5, 1.5)
    maxIter = 100
    no_plot = true
end

optimizer = MOI.instantiate(AB.LazyAbstraction.Optimizer)

Udom = DO.DomainList(Ugrid)
DO.add_set!(Udom, concrete_system.U, DO.OUTER)
AB.LazyAbstraction.set_optimizer!(
    optimizer,
    concrete_problem,
    maxIter,
    pre_image,
    post_image,
    compute_reachable_set,
    minimum_transition_cost,
    hx_heuristic,
    hx,
    Udom;
    periodic_dims = periodic_dims,
    periodic_periods = periodic_periods,
    periodic_start = periodic_start,
)

# Build the state feedback abstraction and solve the optimal control problem using A* algorithm
MOI.optimize!(optimizer)

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
abstract_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_lyap_fun"))
concrete_lyap_fun = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_lyap_fun"))

# ## Simulation
# We define the cost and stopping criteria for a simulation
cost_eval(x, u) = UT.function_value(concrete_problem.transition_cost, x, u)
reached(x) = x ∈ concrete_problem.target_set
nstep = typeof(concrete_problem.time) == PR.Infinity ? 100 : concrete_problem.time; # max num of steps
# We simulate the closed loop trajectory
x0 = UT.get_center(concrete_problem.initial_set)
x_traj, u_traj = ST.get_closed_loop_trajectory(
    concrete_system,
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

c_traj, cost_true = ST.get_cost_trajectory(x_traj, u_traj, cost_eval)
cost_bound = concrete_lyap_fun(x0)
println("Goal set reached")
println("Guaranteed cost:\t $(cost_bound)")
println("True cost:\t\t $(cost_true)")

# ## Display the results
# # Display the specifications and domains
fig1 = plot(; aspect_ratio = :equal)
#We display the concrete domain
plot!(fig1, concrete_system.X; color = :yellow, opacity = 0.5)

#We display the abstract domain
plot!(fig1, abstract_system.Xdom; color = :blue, opacity = 0.5)

#We display the concrete specifications
plot!(fig1, concrete_problem.initial_set; color = :green, opacity = 0.8)
plot!(fig1, concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.8)

#We display the concrete trajectory
plot!(fig1, x_traj; ms = 0.5)

# # Display the abstraction and Lyapunov-like function
fig2 = plot(; aspect_ratio = :equal)
plot!(fig2, abstract_system; dims = [1, 2], cost = true, lyap_fun = optimizer.lyap_fun)

# # Display the Bellman-like value function (heuristic)
fig3 = plot(; aspect_ratio = :equal)
plot!(
    fig3,
    optimizer.abstract_system_heuristic;
    arrowsB = false,
    dims = [1, 2],
    cost = true,
    lyap_fun = optimizer.bell_fun,
)

# # Display the results of the A* algorithm
fig4 = plot(; aspect_ratio = :equal)
plot!(fig4, optimizer.lazy_search_problem)

@testset "lazy_abstraction" begin
    @test cost_bound ≈ 9.0 rtol = 1e-3
    @test cost_true ≈ 9.0 rtol = 1e-3
    @test cost_true <= cost_bound
    @test isa(fig1, Plots.Plot{Plots.GRBackend})
    @test isa(fig2, Plots.Plot{Plots.GRBackend})
    @test isa(fig3, Plots.Plot{Plots.GRBackend})
    @test isa(fig4, Plots.Plot{Plots.GRBackend})
end
end
