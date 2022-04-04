include(joinpath("..", "Abstraction", "abstraction.jl"))
#include(joinpath("..", "CDC2021","CDC2021.jl"))
include(joinpath("..", "CDC2021","general_domain.jl"))
include(joinpath("..", "CDC2021","utils.jl"))
include(joinpath("..", "CDC2021","lazy_abstraction.jl"))
include(joinpath("..", "CDC2021","alternating_simulation.jl"))
module TestM2
using Test, StaticArrays,Plots, LinearAlgebra, LazySets
using LinearAlgebra, Statistics
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

using ForwardDiff, IntervalArithmetic
using HybridSystems
include("sorted_vector_set.jl")
include("automaton.jl")
include("domain.jl")
include("markov_chain.jl")
include("symbolic.jl")
# include("algo.jl")
include("multigrid.jl")
include("system2.jl")


function build_dom()
    hx = [0.5, 0.5, 1.75, 7.0, 7.0] #*0.3
    r1 = 7.0
    r2 = 7.0
    r3 = 7.0
    r4 = 7.0
    r5 = 7.0
    X = AB.HyperRectangle(SVector(-r1, -r2, -r3, -r4, -r5), SVector(r1, r2, r3, r4, r5))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0, -100.0, -100.0, -100.0), SVector(-99.0, -99.0, -99.0, -99.0, -99.0))
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;fit=true)  #il faut fit = true !!!
    return X,Xdom
end
struct SimpleSystem{N,T,F<:Function,F2<:Function} <: AB.ControlSystem{N,T}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F
    f::F2
end
function NewControlSystemGrowthRK4(tstep, F_sys, measnoise::SVector{N,T}, nsys) where {N,T}
    sys_map = let nsys = nsys
        (x::SVector{N,T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
    end
    return SimpleSystem(tstep, measnoise, sys_map, F_sys)
end

function build_system()
    tstep = 0.08
    measnoise = SVector(0.0,0.0,0.0,0.0,0.0)
    sysnoise = SVector(0.0,0.0,0.0,0.0,0.0)
    nsys = 5

    f = HD_f_with_boundary(HD_f(),5.5,-0.8)#-0.8
    contsys = NewControlSystemGrowthRK4(tstep, f, measnoise, nsys)
    return contsys
end

function test_limit_cycle()
    println("BUILD DOMAIN")
    X,Xdom = build_dom()
    Udom = UT.CustomList([SVector(0.0,0.0,0.0,0.0)])
    println("BUILD SYSTEM")
    sys = build_system()
    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    println("COMPUTE TRANSITIONS")
    compute_symbolic_full_domain!(symmodel,sys)
    println("MARKOV CHAIN")
    update_MC!(symmodel)
    println("PLOT")
    plot_steady_state(sys,symmodel,xlims=[-10.0,10.0],ylims=[-10.0,10.0])
    plot_steady_state(sys,symmodel,fact=0.0015)
    plot_steady_state(sys,symmodel,fact=0.001)
    plot_steady_state(sys,symmodel,fact=0.0005)
    plot_steady_state(sys,symmodel,fact=0.0001)
    plot_steady_state(sys,symmodel,fact=0.00001)
    plot_steady_state(sys,symmodel,fact=0.000001)
    plot_steady_state(sys,symmodel,fact=0.0000001)
    plot_steady_state(sys,symmodel,fact=0.00000001)
    plot_steady_state(sys,symmodel,fact=0.000000001)
    plot_steady_state(sys,symmodel,fact=0.0000000001)
    plot_steady_state(sys,symmodel,fact=0.00000000001)
    plot_steady_state(sys,symmodel,fact=0.000000000001)
    plot_steady_state(sys,symmodel,fact=0.0000000000001)
    plot_steady_state(sys,symmodel,fact=0.00000000000001)
    plot_steady_state(sys,symmodel,fact=0.000000000000001)
    plot_steady_state(sys,symmodel,fact=0.0000000000000001)
    plot_steady_state(sys,symmodel,fact=0.00000000000000001)
    plot_steady_state(sys,symmodel,fact=0.000000000000000001)
    plot_steady_state(sys,symmodel,fact=0.0000000000000000001)
    plot_steady_state(sys,symmodel,fact=0.0)
    plot_steady_state(sys,symmodel,fact=0.0005,dims=[1,3])
    plot_steady_state(sys,symmodel,dims=[1,3])
    plot_steady_state(sys,symmodel,dims=[1,3],fact=0.0)
    # plot_steady_state(sys,symmodel,fact=0.0001,xlims=[-5.0,5.0])
    # plot_steady_state(sys,symmodel,fact=0.0,xlims=[-5.0,5.0])
    # plot_shannon_entropy(sys,symmodel,xlims=[-8.0,8.0], ylims=[-5.0,10.0])

    # region = get_steady_state_region_threshold(symmodel, 0.0005)#0.0008
    # plot_region(symmodel,sys,region)
    # fixed_point!(sys,symmodel,region)
    # println("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL")
    # println(region)
    # plot_region(symmodel,sys,region)
end
test_limit_cycle()




function compute_steady_state_Sparse(M)
    x = [1; (I - M[2:end,2:end]) \ Vector(M[2:end,1])]
    x =  x/sum(x)
    a = norm(M * x - x) / norm(x)
    println(a)
    return x
end


function compute_steady_state_Sparse_2(M)
    println(size(M))
    n = size(M)[1]
    L = [I-M;ones(1,n)]
    println(size(L))
    b = [zeros(n,1);1]
    x = L\b
    println(x)
    return x
end

function test()
    # A = sprand(50,50, 1e-2)
    # M = A./ sum(A, dims=2)
    # bool = is_Markov_matrix_Sparse(M)
    # println(bool)
    #
    # M = sparse(M')
    # find an eigenvector x corresponding to eigenvalue λ = 1,
    # normalized so that the first component is 1
    A = sprand(1000,1000, 1e-2)
    A = A ./ sum(A, dims=2)
    compute_steady_state_Sparse(A)
end

function test2()
    I = [1 , 1, 2, 3]
    J = [2, 3 , 1, 1]
    V = [0.5 , 0.5, 1.0, 1.0]
    A = sparse(I,J,V)
    bool = is_Markov_matrix_Sparse(A)
    println(bool)
    compute_steady_state_Sparse_2(A)

end

# test2()
end
