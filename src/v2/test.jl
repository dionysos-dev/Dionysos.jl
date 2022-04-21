include("../Dionysos.jl")
using ..Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const LA = DI.Control.LazyAbstractionReach


using StaticArrays,Plots, DiscreteMarkovChains, LinearAlgebra
using Distributions

using HybridSystems
include("sorted_vector_set.jl")
include("automaton.jl")
include("domain.jl")
include("symbolic.jl")
include("algo.jl")


function build_dom()
    X = AB.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = AB.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [3.0, 1.0]*2.0
    periodic = Int[1]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0,fit=true)
    return X,Xdom
end

function build_dom2()
    X = AB.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = AB.HyperRectangle(SVector(-5.0, -5.0), SVector(-4.0, -4.0)) #AB.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [30.0, 30.0]
    periodic = Int[1,2]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0,fit=true)
    return X,Xdom
end

println("START")

function test1()
    nstates = 4
    nsymbols = 3
    autom = NewProbaAutomaton(nstates, nsymbols)
    println(autom.nstates)
    add_transition!(autom, 1,2,3,0.5)
    add_transition!(autom, 2,1,3,0.1)
    add_transition!(autom, 3,2,3,0.2)
    add_transition!(autom, 4,0,3,0.4)
    add_transition!(autom, 1,2,4,0.4)

    println(ntransitions(autom))
    targetlist= []
    compute_post!(targetlist, autom, 1, 2)
    println(targetlist)

    for e in pre(autom,3)
        println(e)
    end

    delete_transition!(autom, 3,2,3)
    println(ntransitions(autom))
end


function NewSimpleControlSystem(tstep,measnoise::SVector{N,T}) where {N,T}
    function sys_map(x::SVector{N,T}, u, tstep)
        return x+tstep*u
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end

function build_system()
    tstep = 7.5#0.8
    measnoise = SVector(0.0, 0.0)
    return NewSimpleControlSystem(tstep,measnoise)
end

function NewSimpleControlSystem2(tstep,measnoise::SVector{N,T}) where {N,T}
    function sys_map(x::SVector{N,T}, u, tstep)
        θ = π/6
        x1 = cos(θ)*(x[1]-15.0) - sin(θ)*(x[2]-15.0)
        x2 = sin(θ)*(x[1]-15.0) + cos(θ)*(x[2]-15.0)
        y = SVector(x1+15.0,x2+15.0)
        return y
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end

function build_system2()
    tstep = 0.8 #0.8
    measnoise = SVector(0.0, 0.0)
    return NewSimpleControlSystem2(tstep,measnoise)
end

function post_image(symmodel,sys,xpos,l,u)
    dom = symmodel.Xdom.domains[l]
    x = AB.get_coord_by_pos(dom.grid, xpos)
    tstep = sys.tstep
    Fx = sys.sys_map(x, u, tstep)
    r = dom.grid.h/2.0 + sys.measnoise
    Fr = r
    rectI = AB.get_pos_lims_outer(dom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(AB._ranges(rectI)...)
    over_approx = Int[]
    allin = true #but assume domain is forward invariant mais necessaire +- pour les obstacles...
    for ypos in ypos_iter
        ypos = D.set_in_period_pos(dom,ypos)
        if !(ypos in dom)
            allin = false
            break
        end
        target = get_state_by_xpos(symmodel, ypos,l)
        push!(over_approx, target)
    end
    sort!(over_approx)
    unique!(over_approx)
    return allin ? over_approx : []
end

function test6()
     P = [1.0 0.0 0.0; 0.0 1.0 0.0;0.0 0.0 1.0]
     transition_matrix = [
    0.0 1.0 0.0;
    1.0 0.0 0.0;
    0.0 0.0 1.0
    ]
    chain = DiscreteMarkovChain(transition_matrix)
    p = stationary_distribution(chain)
    println(p)

end

function test7()
    rec = AB.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 5.0))

    #sample(rec,4)
    tab = [5,9,8,6,5,4,5,14,12,14,15]
    println(sort(tab))
    r = count_occurences(tab)
    println(r)
end

function NewSimpleControlSystem3(tstep,measnoise::SVector{N,T}) where {N,T}
    function sys_map(x::SVector{N,T}, u, tstep)
        c = SVector(15,15)
        D = @SMatrix [ 1.0 0.0 ;
                       0.0  1.0 ]

        a = (0.9-0.5)/20.0
        f = 0.95


        if x[1]>15.0 && x[2]<15.0
            f = 1.1
        end

        if norm(x-c,2)>12.0
            f = 0.6
        end
        D = D*f
        B = @SMatrix [1.0 0.0 ;
                      0.0 1.0]
        C = inv(B)
        θ = π/6
        E = @SMatrix [ cos(θ) -sin(θ) ;
                       sin(θ)  cos(θ)]
        A = D*E
        # B = Array(A)
        # eigs = eigvals(B)
        # println(abs(eigs[1]))
        # println(abs(eigs[2]))
        # println(eigs)
        # println(A)
        # println(A)
        # println(x)

        #println(x-c)
        y = A*(x-c) + c
        #println(y)
        return y
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end

function build_system3()
    tstep = 0.8 #0.8
    measnoise = SVector(0.0, 0.0)
    return NewSimpleControlSystem3(tstep,measnoise)
end



function RungeKutta4(F, x, u, tstep, nsub::Int)
    τ = tstep/nsub
    for i in 1:nsub
        Fx1 = F(x, u)
        xrk = x + Fx1*(τ/2.0)
        Fx2 = F(xrk, u)
        xrk = x + Fx2*(τ/2.0)
        Fx3 = F(xrk, u)
        xrk = x + Fx3*τ
        Fx4 = F(xrk, u)
        x = x + (Fx1 + Fx2*2.0 + Fx3*2.0 + Fx4)*(τ/6.0)
    end
    return x
end

function NewControlSystemGrowthRK4(tstep, F_sys, measnoise::SVector{N,T}, nsys) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end


function build_system5()
    tstep = 0.08#0.8
    measnoise = SVector(0.0, 0.0)
    nsys = 4
    function F_sys(x::SVector{N,T}, u) where {N,T}
        μ = 0.0
        y1 = x[2]
        y2 = μ*(1-x[1]*x[1])*x[2]-x[1]
        y = SVector(y1,y2)
        # D = @SMatrix [ 1.0 0.0 ;
        #                0.0  1.0 ]
        # A = @SMatrix [ 0.0 -1.0 ;
        #              1.0 0.0 ]
        # a = -0.08
        # b = 1.0
        # if x[1]>0.0 && x[2]<0.0
        #     a = 1.5
        #     b = 3.0
        # end
        if norm(x,2)>8.0
            a = -0.8
            b = 0.0#1.0
            A = @SMatrix [ a -b ;
                           b a ]
            y = A*x
        end
        return y
    end
    return NewControlSystemGrowthRK4(tstep,F_sys,measnoise,nsys)
end

function build_system4()
    tstep = 0.8#0.8
    measnoise = SVector(0.0, 0.0)
    nsys = 4
    function F_sys(x::SVector{N,T}, u) where {N,T}
        D = @SMatrix [ 1.0 0.0 ;
                       0.0  1.0 ]
        A = @SMatrix [ 0.0 -1.0 ;
                     1.0 0.0 ]
        a = -0.08
        b = 1.0
        if x[1]>0.0 && x[2]<0.0
            a = 1.5
            b = 3.0
        end
        if norm(x,2)>40.0
            a = -0.5
            b = 1.0
        end
        A = @SMatrix [ a -b ;
                       b a ]

        y = A*x
        #println(y)
        return y
    end
    return NewControlSystemGrowthRK4(tstep,F_sys,measnoise,nsys)
end

function build_dom4()
    X = AB.HyperRectangle(SVector(-50.0, -50.0), SVector(50.0, 50.0))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0)) #AB.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [5.0, 5.0]#*2.0
    periodic = Int[]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0,fit=true)
    return X,Xdom
end

function build_dom5()
    X = AB.HyperRectangle(SVector(-10.0, -10.0), SVector(10.0, 10.0))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0)) #AB.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [1.0, 1.0]*0.5#*2.0
    periodic = Int[]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0,fit=true)
    return X,Xdom
end

function test10()
    X,Xdom = build_dom4()
    Udom = UT.CustomList([SVector(0.0,1.0)]) #1.0,1.0
    sys = build_system4()


    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    obs = AB.HyperRectangle(SVector(7.0, 0.0), SVector(14.0, 10.0))#AB.HyperRectangle(SVector(7.0, 0.0), SVector(14.0, 10.0))
    prob = NewSafetyProblem(symmodel,sys,post_image,obs,false)
    solve!(prob)
    println("PLOTTING")
    #fig = plot(aspect_ratio = 1,legend = false)
    #plot_automaton!(prob)
    # for i in 1:get_ncells(prob.symmodel)
    #     if prob.active[i]
    #         plot_automaton(prob,bool=true,src=i)
    #     end
    # end

    #plot_automaton(prob)
end

function build_dom3()
    X = AB.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = AB.HyperRectangle(SVector(-5.0, -5.0), SVector(-4.0, -4.0)) #AB.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [5.0, 5.0]*6.0#*2.0
    periodic = Int[]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0,fit=true)
    return X,Xdom
end

function test8()
    sys = build_system5()
    x0 = SVector(1.0,1.0) #SVector(20.0,20.0)
    l = 100

    center = [0.0,0.0]
    r = [50.0,50.0]
    fig = plot(aspect_ratio = 1,legend = false)
    #plot!(D.rectangle(center,r), opacity=.9,color=:yellow)
    plot_trajectory!(sys,x0,l)
    display(fig)
end
function test9()
    center = [0.0,0.0]# [15.0,15.0]
    r = [50.0,50.0]#[15.0,15.0]
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(D.rectangle(center,r), opacity=.9,color=:yellow)
    sys = build_system4()
    l = 30
    N = 20
    for k in 1:N
        x = k*50.0/N
        y = x
        tab = [x,y]
        x0 =  SVector{2,Float64}(tab)
        plot_trajectory!(sys,x0,l)
    end
    display(fig)
end
function test15()
    center = [0.0,0.0]# [15.0,15.0]
    r = [10.0,10.0]#[15.0,15.0]
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(D.rectangle(center,r), opacity=.9,color=:yellow)
    sys = build_system5()
    l = 100
    N = 5
    for k in 1:N
        x = k*5.0/N
        y = x
        tab = [x,y]
        x0 =  SVector{2,Float64}(tab)
        plot_trajectory!(sys,x0,l)
    end
    display(fig)
end

function test5()
    X,Xdom = build_dom5()
    Udom = UT.CustomList([SVector(0.0,1.0)]) #1.0,1.0
    sys = build_system5()


    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    obs = AB.HyperRectangle(SVector(7.0, 0.0), SVector(14.0, 10.0))#AB.HyperRectangle(SVector(7.0, 0.0), SVector(14.0, 10.0))
    prob = NewSafetyProblem(symmodel,sys,post_image,obs,false)
    solve!(prob)
    println("PLOTTING")
    #fig = plot(aspect_ratio = 1,legend = false)
    #plot_automaton!(prob)
    # for i in 1:get_ncells(prob.symmodel)
    #     if prob.active[i]
    #         plot_automaton(prob,bool=true,src=i)
    #     end
    # end

    #plot_automaton(prob)
end


# test8()
test5()

end # end module
