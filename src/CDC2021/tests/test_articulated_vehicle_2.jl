include(joinpath("..", "..", "Abstraction", "abstraction.jl"))
include("../general_domain.jl")
include("../utils.jl")
include("../alternating_simulation.jl")
include("../partition.jl")
include("../lazy_abstraction.jl")
include("../branch_and_bound.jl")
include("../optimal_control.jl")

module TestArticuledVehicle
using ForwardDiff, IntervalArithmetic
using LinearAlgebra,StaticArrays,Random
using Plots, StaticArrays, JuMP

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const U = Utils

using ..BranchAndBound
const BB = BranchAndBound

using ..OptimalControl
const OC = OptimalControl

using ..Lazy_abstraction
const LA = Lazy_abstraction

using ..AlternatingSimulation
const AS = AlternatingSimulation

const α = 8.0

function plot_vehicle!(x,u)
    function plot_wheel!(p,θ,l,r)
        y1 = [p[1] + l*cos(θ-π/2.0), p[2] + l*sin(θ-π/2.0)]
        y2 = [p[1] - l*cos(θ-π/2.0), p[2] - l*sin(θ-π/2.0)]
        plot!([y1[1],y2[1]], [y1[2],y2[2]],color =:black,linewidth = 4)
        y3 = [y1[1] + r*cos(θ), y1[2] + r*sin(θ)]
        y4 = [y1[1] - r*cos(θ), y1[2] - r*sin(θ)]
        plot!([y3[1],y4[1]], [y3[2],y4[2]],color =:black,linewidth = 4)
        y5 = [y2[1] + r*cos(θ), y2[2] + r*sin(θ)]
        y6 = [y2[1] - r*cos(θ), y2[2] - r*sin(θ)]
        plot!([y5[1],y6[1]], [y5[2],y6[2]],color =:black,linewidth = 4)
    end
    L1 = 1.0*α
    L2 = 0.5*α #1.5 ??
    c = 0.3*α
    l = 0.25*α
    r = 0.15*α
    x1 = [x[1] + L1*cos(x[3]), x[2] + L1*sin(x[3])]
    x2 = [x[1] - c*cos(x[3]),  x[2] - c*sin(x[3])]
    x3 = [x2[1] - L2*cos(x[3]+x[4]),  x2[2] - L2*sin(x[3]+x[4])]
    plot!([x1[1],x2[1]], [x1[2],x2[2]],color =:black,linewidth = 4)
    plot!([x3[1],x2[1]], [x3[2],x2[2]],color =:black,linewidth = 4)
    scatter!([x2[1]],[x2[2]],color =:red,markersize=4)
    v,δ = u
    x4 = [x1[1] + r*cos(x[3]+δ), x1[2] + r*sin(x[3]+δ)]
    x5 = [x1[1] - r*cos(x[3]+δ), x1[2] - r*sin(x[3]+δ)]
    plot!([x4[1],x5[1]], [x4[2],x5[2]],color =:black,linewidth = 4)
    plot_wheel!([x[1],x[2]],x[3],l,r)
    plot_wheel!(x3,x[3]+x[4],l,r)
end

function plot_articulted_vehicle_traj(traj;Xdom=nothing)
    for e in traj
        x,u = e
        fig = plot(aspect_ratio = 1,legend = false,xlims = (-4,4).*α,ylims = (-4,4).*α)
        if Xdom != nothing
            U.plot_domain!(Xdom,dims=[1,2])
        end
        plot_vehicle!(x,u)
        sleep(0.02)
        display(fig)
    end
end
function test3()
    tstep = 0.1
    X,O = build_domain()
    contsys = build_system(X)
    x = SVector(0.0,0.0,π/2.0,0.0)
    u = [1.0,π/3.0]

    fig = plot(aspect_ratio = 1,legend = false,xlims = (0,10))
    traj = []
    for i=1:60
        push!(traj,(x,u))
        fig = plot(aspect_ratio = 1,legend = false,xlims = (-3,2).*α,ylims = (-2,2).*α)
        plot_vehicle!(x,u)
        x = contsys.sys_map(x, u, tstep)
        sleep(0.02)
        display(fig)
    end
    plot_articulted_vehicle_traj(traj)
end

## define the system

function NewControlSystemGrowthLx(tstep, f, sysnoise::SVector{N,T}, measnoise::SVector{N,T}, nsys, ngrowthbound,X) where {N,T}
    # Use `let` following the advice of
    # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured

    function L_growthbound(u,K)
        #L = ForwardDiff.jacobian(x->f(x,u),K)
        L = jacobian(u,K)
        M = @SMatrix [i == j ? max(L[i,j].lo,L[i,j].hi) : max(abs(L[i,j].lo), abs(L[i,j].hi)) for i in 1:4, j in 1:4]
        return M
    end
    sys_map = let nsys = nsys
        (x, u, tstep) ->
            AB.RungeKutta4(f, x, u, tstep, nsys)
    end
    function growthbound_map(r::SVector{N,T}, u, tstep, x; nstep=1)
        K = compute_K(X,r,u,tstep,x, nstep=nstep)
        L = L_growthbound(u, K)
        function F_growthbound(r, u)
            return L*r + sysnoise
        end
        return AB.RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)
    end
    return AB.ControlSystemGrowth(tstep, sysnoise, measnoise, sys_map, growthbound_map)
end

function jacobian(u,x)
    L1 = 1.0*α
    L2 = 0.5*α
    c = 0.25*α
    v, δ = u
    tanδ = tan(δ)
    a = -(v/(L1*L2))*(L1*cos(x[4])-c*sin(x[4])*tanδ)
    return SMatrix{4,4}(0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, -v*sin(x[3]),v*cos(x[3]),0.0,0.0, 0.0,0.0,0.0,a)
end

function compute_K(X,r,u,tstep,x;nstep=1)
    K = IntervalBox(X.lb, X.ub)
    R = IntervalBox(-r,r)
    for k=1:nstep
        B = f(K,u)
        K = x + R + B*tstep
        K = K ∪ (x+R)
    end
    return K
end

function f(x, u)
    v, δ = u
    tanδ = tan(δ)
    L1 = 1.0*α
    L2 = 0.5*α
    c = 0.3*α
    sθ, cθ = sincos(x[3])
    dx = v * cθ
    dy = v * sθ
    dθ = v * tanδ / L1
    sφ, cφ = sincos(x[4])
    dφ = -v * (L1*sφ + (c * cφ + L2) * tanδ) / (L1 * L2)
    return SVector(dx, dy, dθ, dφ)
end

## problem specific-functions required
function minimum_transition_cost(symmodel,contsys,source,target)
    return 1.0
end

function compute_reachable_set(rect::AB.HyperRectangle{SVector{N,T}},contsys,Udom) where {N,T}
    tstep = contsys.tstep
    r = (rect.ub-rect.lb)/2.0 + contsys.measnoise
    x = U.center(rect)
    n =  U.dims(rect)
    lb = SVector(ntuple(i -> Inf, Val(N)))
    ub = SVector(ntuple(i -> -Inf, Val(N)))
    for upos in AB.enum_pos(Udom)
        u = AB.get_coord_by_pos(Udom.grid,upos)
        Fx = contsys.sys_map(x, u, tstep)
        Fr = contsys.growthbound_map(r, u, tstep, x, nstep=5)
        lb = min.(lb,Fx .- Fr)::SVector{N,T}
        ub = max.(ub,Fx .+ Fr)::SVector{N,T}
    end
    return AB.HyperRectangle(lb,ub)
end


function post_image(symmodel,contsys,xpos,u)
    Xdom = symmodel.Xdom
    x = AB.get_coord_by_pos(Xdom.grid, xpos)
    tstep = contsys.tstep
    Fx = contsys.sys_map(x, u, tstep)
    r = Xdom.grid.h/2.0 + contsys.measnoise
    Fr = contsys.growthbound_map(r, u, tstep, x, nstep=5)
    post_rec = AB.HyperRectangle(Fx .- Fr, Fx .+ Fr)
    rectI = AB.get_pos_lims_outer(Xdom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(AB._ranges(rectI)...)
    over_approx = Int[]
    for ypos in ypos_iter
        ypos = D.set_in_period_pos(Xdom,ypos)
        if !(ypos in Xdom)
            empty!(over_approx)
            break
        end
        target = AB.get_state_by_xpos(symmodel, ypos)
        push!(over_approx, target)
    end
    return over_approx
end


function pre_image(symmodel,contsys,xpos,u)
    return post_image(symmodel,contsys,xpos,SVector(-u[1],u[2]))
end

##
function build_domain()
    X = AB.HyperRectangle(SVector(0.0,0.0,-π,-π/3.0),SVector(60.0,60.0,π,π/3.0))
    O1 = AB.HyperRectangle(SVector(40.1,0.0,-π,-π),SVector(60.0,19.9,π,π))
    O2 = AB.HyperRectangle(SVector(40.1,40.1,-π,-π),SVector(60.0,60.0,π,π))
    O = [O1,O2]
    return X,O
end

function build_system(X)
    tstep = 1.8
    sysnoise = SVector(0.0, 0.0, 0.0, 0.0)
    measnoise = SVector(0.0, 0.0, 0.0, 0.0)
    nsys = 5
    ngrowthbound = 5
    contsys = NewControlSystemGrowthLx(tstep, f, sysnoise, measnoise, nsys, ngrowthbound, X)
    return contsys
end


function build_dom()
    X,obstacle = build_domain()
    hx = [1.0, 1.0, 0.1, 0.1]
    periodic = [3,4]
    periods = [2*π,2*π]
    T0 = [-π,-π]
    grid = D.build_grid_in_rec(X, hx)
    d = D.RectanglularObstacles(X, obstacle)
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0)
    return X,Xdom
end

function build_Udom()
    U = AB.HyperRectangle(SVector(-1.0,-40/180*π), SVector(1.0,40/180*π))
    x0 = SVector(-1.0,0.0); hu = SVector(2.0,(10.0/180.0)*π)
    Ugrid = AB.GridFree(x0,hu)
    Udom = AB.DomainList(Ugrid)
    AB.add_set!(Udom, U, AB.OUTER)
    println(AB.get_ncells(Udom))
    #=for pos in AB.enum_pos(Udom)
        u = AB.get_coord_by_pos(Udom.grid,pos)
        println(u)
    end
    println(AB.get_ncells(Udom))
    println(Udom)=#
    return Udom
end
function transition_cost(x,u)
    return 1.0
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
    hx = [5.0, 5.0, 0.4, 0.4]
    periodic = [3,4]
    periods = [2*π,2*π]
    T0 = [-π,-π]
    Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods,T0=T0)
    AB.add_set!(Xdom, X , AB.OUTER)
    symmodel = AB.NewSymbolicModelListList(Xdom, Udom)
    problem = AS.symmodelProblem(symmodel,contsys,compute_reachable_set,minimum_transition_cost,AS.get_possible_transitions_2)
    autom = AS.build_alternating_simulation(problem)
    symmodel = AB.with_automaton(symmodel, autom)
    # build the heurisitic
    initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
    heuristic_data = AS.build_heuristic(symmodel,initlist)
    println("Heuristic ended")
    println("kkk")
    #fig = plot(aspect_ratio = 1,legend = false)
    #U.plot_domain!(heuristic_data.symmodel.Xdom)
    #AS.plot_heuristic!(heuristic_data)
    #display(fig)
    return heuristic_data
end

function test2()
    println()

    contsys = build_system(X)
    Udom = build_Udom()
    periodic = [3,4]
    periods = [2*π,2*π]
    T0 = [-π,-π]
    x = SVector(0.01,0.01,π/2.0+0.2,0.0)
    tstep = 1.8
    u = [-1.0,80/180*π]#[3.0,π/3.0]
    h = SVector(20.0, 20.0, 1.2, 2*π)
    h = SVector(1.0, 1.0, 0.05, 0.05)
    r = h/2
    rec = AB.HyperRectangle(x .- r, x .+ r)
    # post-image
    Fx = contsys.sys_map(x, u, tstep)
    Fr = contsys.growthbound_map(r, u, tstep, x; nstep=5)
    R = AB.HyperRectangle(Fx .- Fr, Fx .+ Fr)
    #post-K
    K = compute_K(X,r,u,tstep,x,nstep=6)
    n = 4; lb = zeros(n); ub = zeros(n)
    for i=1:n
        lb[i] = K[i].lo
        ub[i] = K[i].hi
    end
    K = AB.HyperRectangle(lb,ub)

    #pre-K
    K2 = compute_K(X,r,u,-tstep,x,nstep=5)
    n = 4; lb = zeros(n); ub = zeros(n)
    for i=1:n
        lb[i] = K2[i].lo
        ub[i] = K2[i].hi
    end
    K2 = AB.HyperRectangle(lb,ub)
    #reachable set
    RS = compute_reachable_set(rec,contsys,Udom)
    # pre-image
    u2 = SVector(-u[1],u[2])
    Fx2 = contsys.sys_map(x, u2, tstep)
    Fr2 = contsys.growthbound_map(r, u2, tstep, x; nstep=5)
    R2 = AB.HyperRectangle(Fx2 .- Fr2, Fx2 .+ Fr2)

    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(RS.lb[1:2],RS.ub[1:2]),color=:black)
    plot!(U.rectangle(K.lb[1:2],K.ub[1:2]),color=:green)
    plot!(U.rectangle(K2.lb[1:2],K2.ub[1:2]),color=:yellow)
    plot!(U.rectangle(rec.lb[1:2],rec.ub[1:2]),color=:blue)
    plot!(U.rectangle(R.lb[1:2],R.ub[1:2]),color=:red)
    plot!(U.rectangle(R2.lb[1:2],R2.ub[1:2]),color=:pink)
    display(fig)
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(RS.lb[3:4],RS.ub[3:4]),color=:black)
    plot!(U.rectangle(K.lb[3:4],K.ub[3:4]),color=:green)
    plot!(U.rectangle(K2.lb[3:4],K2.ub[3:4]),color=:yellow)
    plot!(U.rectangle(rec.lb[3:4],rec.ub[3:4]),color=:blue)
    plot!(U.rectangle(R.lb[3:4],R.ub[3:4]),color=:red)
    plot!(U.rectangle(R2.lb[3:4],R2.ub[3:4]),color=:pink)
    display(fig)
end



function test()

    println("start")
    X,O = build_domain()
    contsys = build_system(X)

    Udom = build_Udom()


    X,Xdom = build_dom()
    Udom = build_Udom()
    # build system
    symmodel = D._SymbolicModel(Xdom, Udom)

    # control problem
    x0 = SVector(30.0, 5.0,π/2.0,0.0)
    RI = SVector(0.5, 0.5, 0.005, 0.005)
    _I_ = AB.HyperRectangle(SVector(30.0, 5.0,π/2.0,0.0)-RI, SVector(30.0, 5.0,π/2.0,0.0)+RI)
    #_T_ = AB.HyperRectangle(SVector(25.0, 12.0,π/2.0-0.1,-0.1), SVector(35.0, 19.0,π/2.0+0.1,0.1))
    _T_ = AB.HyperRectangle(SVector(25.0, 12.0,π/2.0-0.3,-0.3), SVector(35.0, 19.0,π/2.0+0.3,0.3))
    initlist = U.get_symbol(symmodel,_I_,AB.OUTER)
    targetlist = U.get_symbol(symmodel,_T_,AB.INNER)
    # Heuristic data
    println(initlist)
    println(targetlist)
    heuristic_data = build_heuristic_data(X,contsys,Udom,_I_)
    println()
    println()
    println("lll")
    println()
    println()
    println(initlist)
    println("LLLLM%%%%%%%%%%%")
    println(targetlist)
    # Lazy Abstraction implementation
    time = @elapsed begin
    problem,sucess = LA.compute_controller(symmodel, contsys, initlist, targetlist, transition_cost, pre_image, post_image, h, heuristic_data=heuristic_data)
    contr = problem.contr
    end
    println("total time: lazy abstraction + controller: ", time)
    fig = plot(aspect_ratio = 1,legend = false)
    LA.plot_result!(problem,x0=x0)
    display(fig)
    fig = plot(aspect_ratio = 1,legend = false)
    LA.plot_result!(problem,dims=[3,4],x0=x0)
    display(fig)
end


test()

end # end module
