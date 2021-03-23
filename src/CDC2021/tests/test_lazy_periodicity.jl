include(joinpath("..", "..", "Abstraction", "abstraction.jl"))
include("../general_domain.jl")
include("../utils.jl")
include("../lazy_abstraction.jl")
include("../alternating_simulation.jl")

module TestMain
using Test, StaticArrays,Plots

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const UT = Utils

using ..Lazy_abstraction
const LA = Lazy_abstraction

using ..AlternatingSimulation
const AS = AlternatingSimulation



## specific function
function post_image(symmodel,contsys,xpos,u)
    Xdom = symmodel.Xdom
    x = AB.get_coord_by_pos(Xdom.grid, xpos)
    tstep = contsys.tstep
    Fx = contsys.sys_map(x, u, tstep)
    r = Xdom.grid.h/2.0 + contsys.measnoise
    Fr = r

    rectI = AB.get_pos_lims_outer(Xdom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(AB._ranges(rectI)...)
    over_approx = []
    allin = true
    for ypos in ypos_iter
        ypos = D.set_in_period_pos(Xdom,ypos)
        if !(ypos in Xdom)
            allin = false
            break
        end
        target = AB.get_state_by_xpos(symmodel, ypos)
        push!(over_approx, target)
    end
    return allin ? over_approx : []
end

function pre_image(symmodel,contsys,xpos,u)
    grid = symmodel.Xdom.grid
    x = AB.get_coord_by_pos(grid, xpos)
    tstep = contsys.tstep
    potential = Int[]
    x_prev = x-tstep*u
    xpos_cell = AB.get_pos_by_coord(grid,x_prev)
    n = 2
    for i=-n:n
        for j=-n:n
            x_n = (xpos_cell[1]+i,xpos_cell[2]+j)
            x_n = D.set_in_period_pos(symmodel.Xdom,x_n)
            if x_n in symmodel.Xdom
                cell = AB.get_state_by_xpos(symmodel, x_n)
                if !(cell in potential)
                    push!(potential,cell)
                end
            end
        end
    end
    return potential
end

function compute_reachable_set(rect::AB.HyperRectangle,contsys,Udom)
    tstep = contsys.tstep
    r = (rect.ub-rect.lb)/2.0 + contsys.measnoise
    Fr = r
    x = UT.center(rect)
    n =  UT.dims(rect)
    lb = fill(Inf,n)
    ub = fill(-Inf,n)
    for upos in AB.enum_pos(Udom)
        u = AB.get_coord_by_pos(Udom.grid,upos)
        Fx = contsys.sys_map(x, u, tstep)
        lb = min.(lb,Fx .- Fr)
        ub = max.(ub,Fx .+ Fr)
    end
    lb = SVector{n}(lb)
    ub = SVector{n}(ub)
    return AB.HyperRectangle(lb,ub)
end

function minimum_transition_cost(symmodel,contsys,source,target)
    return 1.0
end

## poblem data
struct SimpleSystem{N,T,F<:Function} <: AB.ControlSystem{N,T}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F
end

function NewSimpleControlSystem(tstep,measnoise::SVector{N,T}) where {N,T}
    function sys_map(x::SVector{N,T}, u, tstep)
        return x+tstep*u
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end

function build_system()
    tstep = 0.8
    measnoise = SVector(0.0, 0.0)
    return NewSimpleControlSystem(tstep,measnoise)
end

function build_dom()
    X = AB.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = AB.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [0.5, 0.5]
    periodic = Int[]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    grid = D.build_grid_in_rec(X, hx)
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0)
    fig = plot(aspect_ratio = 1,legend = false)
    return X,Xdom
end

function build_Udom()
    U = AB.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
    x0 = SVector(0.0, 0.0)
    hu = SVector(0.5,0.5) #0.5
    Ugrid = AB.GridFree(x0,hu)
    Udom = AB.DomainList(Ugrid)
    AB.add_set!(Udom, U, AB.OUTER)
    box = AB.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))
    AB.remove_set!(Udom, box, AB.OUTER)
    return Udom
end

function transition_cost(x,u)
    return 0.5#1.0
end

## Heursitic
function h(node::LA.S.Node,problem::LA.LazyAbstraction)
    source = node.state.source
    symmodel = problem.symmodel
    xpos = AB.get_xpos_by_state(symmodel, source)
    x = AB.get_coord_by_pos(symmodel.Xdom.grid, xpos)

    heuristic = problem.heuristic_data
    symmodel2 = heuristic.symmodel
    xpos2 = AB.get_pos_by_coord(symmodel2.Xdom.grid, x)
    source2 = AB.get_state_by_xpos(symmodel2, xpos2)
    return heuristic.dists[source2]
end

function build_heuristic_data(X,contsys,Udom,_I_)
    # build the alternating simulation
    hx = [1.0, 1.0]*1.5
    periodic = Int[]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods,T0=T0)
    AB.add_set!(Xdom, X , AB.OUTER)
    symmodel = AB.NewSymbolicModelListList(Xdom, Udom)
    problem = AS.symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,AS.get_possible_transitions_2)
    autom = AS.build_alternating_simulation(problem)
    symmodel = AB.with_automaton(symmodel, autom)
    # build the heurisitic
    initlist = UT.get_symbol(symmodel,_I_,AB.OUTER)
    heuristic_data = AS.build_heuristic(symmodel,initlist)
    println("Heuristic ended")
    #fig = plot(aspect_ratio = 1,legend = false)
    UT.plot_domain!(heuristic_data.symmodel.Xdom)
    AS.plot_heuristic!(heuristic_data)
    #display(fig)
    return heuristic_data
end
##
function test()
    X,Xdom = build_dom()
    Udom = build_Udom()
    # build system
    contsys = build_system()
    symmodel = D._SymbolicModel(Xdom, Udom)
    # control problem
    _I_ = AB.HyperRectangle(SVector(5.0, 5.0), SVector(6.0, 6.0))
    initlist = UT.get_symbol(symmodel,_I_,AB.OUTER)
    _T_ = AB.HyperRectangle(SVector(25.0, 25.0), SVector(28.0, 28.0))
    targetlist = UT.get_symbol(symmodel,_T_,AB.INNER)
    # # Heuristic data
    fig = plot(aspect_ratio = 1,legend = false)
    heuristic_data = build_heuristic_data(X,contsys,Udom,_I_)

    # Lazy Abstraction implementation
    time = @elapsed begin
    problem,sucess = LA.compute_controller(symmodel, contsys, initlist, targetlist, transition_cost, pre_image, post_image, h, heuristic_data=heuristic_data)
    contr = problem.contr
    end
    println("total time: lazy abstraction + controller: ", time)
    x0 = SVector(5.5,5.5)
    LA.plot_result!(problem,x0=x0)
    display(fig)
end
println()
test()
end # end module
