include(joinpath("..", "Abstraction", "abstraction.jl"))
include(joinpath("..", "CDC2021","general_domain.jl"))
include(joinpath("..", "CDC2021","utils.jl"))

module TestM2
using Test, StaticArrays,Plots, LinearAlgebra
using Distributions

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const UT = Utils

using HybridSystems
include("sorted_vector_set.jl")
include("automaton.jl")
include("domain.jl")
include("markov_chain.jl")
include("symbolic.jl")
include("multigrid.jl")

struct SimpleSystem{N,T,F<:Function}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F
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

function NewSimpleControlSystem(tstep,measnoise::SVector{N,T}) where {N,T}
    function sys_map(x::SVector{N,T}, u, tstep)
        return x+tstep*u
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end

function build_system()
    tstep = 7.5
    measnoise = SVector(0.0, 0.0)
    return NewSimpleControlSystem(tstep,measnoise)
end

function build_dom()
    X = AB.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = AB.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [1.0, 1.0]*10.0
    periodic = Int[2]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;periodic=periodic,periods=periods,T0=T0,fit=true)
    return X,Xdom
end

function NewControlSystemGrowthRK4(tstep, F_sys, measnoise::SVector{N,T}, nsys) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end

function test()
    X,Xdom = build_dom()
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys = build_system()
    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    compute_symbolic_full_domain!(symmodel,sys)
    plot_automaton(symmodel)
    update_MC!(symmodel)

    split_node!(symmodel,sys,2)
    update_MC!(symmodel)
    #symmodel.mc = build_Markov_Chain2(symmodel)
    plot_automaton(symmodel)
end

test()

end
