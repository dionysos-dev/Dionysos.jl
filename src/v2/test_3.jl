include(joinpath("..", "Abstraction", "abstraction.jl"))
#include(joinpath("..", "CDC2021","CDC2021.jl"))
include(joinpath("..", "CDC2021","general_domain.jl"))
include(joinpath("..", "CDC2021","utils.jl"))
include(joinpath("..", "CDC2021","lazy_abstraction.jl"))
include(joinpath("..", "CDC2021","alternating_simulation.jl"))


# test lazy abstraction with the simple 2D example (should be deleted as particular case of periodic (general) domain)
module TestM2
using Test, StaticArrays,Plots, DiscreteMarkovChains, LinearAlgebra, LazySets
using Distributions

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


using HybridSystems
include("sorted_vector_set.jl")
include("automaton.jl")
include("domain.jl")
include("markov_chain.jl")
include("symbolic.jl")
include("algo.jl")
include("system.jl")



struct System
    f
    Jacobian
end

function build_dom(f,fi,hx,r)
    X = AB.HyperRectangle(SVector(-r, -r), SVector(r, r))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0))
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d,f=f,fi=fi;fit=true)
    return X,Xdom
end

function plot_norm2_Jacobian(symmodel,Jacobian)
    function norm2Jacobian(x,u)
        return norm(Jacobian(x,u),2)
    end
    plot_map(symmodel,norm2Jacobian)
end



function test_map()
    f = nothing
    fi = nothing
    hx = [1.0, 1.0]*0.3
    X,Xdom = build_dom(f,fi,hx,10.0)
    Udom = UT.CustomList([SVector(0.0,1.0)])
    system = Van_der_pol_oscillator(2.0)
    system_bounded = system_with_boundary(system,8.0,-0.8)
    sys = NewControlSystemGrowthRK4(0.08, system_bounded.f)


    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    prob = NewAnalysisProblem(symmodel,sys,nothing,false)
    compute_symbolic_full_domain!(prob,symmodel,sys)
    update_MC!(prob.symmodel)
    println("Plotting")
    # plot_vector_field(symmodel,system_bounded.f)
    # plot_norm2_Jacobian(symmodel,system_bounded.Jacobian)
    plot_Jacobian(sys,symmodel,system_bounded.Jacobian)
end


test_map()

end
