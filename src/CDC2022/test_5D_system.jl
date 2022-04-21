module Test
include("../Dionysos.jl")
using Test, StaticArrays, Plots, LinearAlgebra, Statistics, Distributions
using ..Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const LA = DI.Control.LazyAbstractionReach

using Test, StaticArrays,Plots, LinearAlgebra
using Distributions

include("nested_domain.jl")
include("proba_automaton.jl")
include("monteCarlo.jl")
include("markov_chain.jl")
include("nested_symbolic.jl")
include("system5D.jl")

function build_dom()
    hx = [0.5, 0.5, 1.75, 7.0, 7.0]
    r1 = 7.0
    r2 = 7.0
    r3 = 7.0
    r4 = 7.0
    r5 = 7.0
    X = UT.HyperRectangle(SVector(-r1, -r2, -r3, -r4, -r5), SVector(r1, r2, r3, r4, r5))
    obstacle = UT.HyperRectangle(SVector(-100.0, -100.0, -100.0, -100.0, -100.0), SVector(-99.0, -99.0, -99.0, -99.0, -99.0))
    d = DO.RectanglularObstacles(X, [obstacle])
    Xdom = DO.GeneralDomainList(hx;elems=d,fit=true)  #il faut fit = true !!!
    return X,Xdom
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
    Udom = DO.CustomList([SVector(0.0,0.0,0.0,0.0)])
    println("BUILD SYSTEM")
    sys = build_system()
    param = Param(5)
    symmodel = NewNestedSymbolicModel(Xdom, Udom, param)
    println("COMPUTE TRANSITIONS")
    compute_symbolic_full_domain!(symmodel,sys)
    println("MARKOV CHAIN")
    update_MC!(symmodel)
    println("PLOT")
    x0 = SVector(2.0,2.0,4.0,4.0,4.0)
    plot_steady_state(sys,symmodel,xlims=[-10.0,10.0],ylims=[-10.0,10.0],x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0015,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0005,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.00001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.00000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.00000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.00000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.000000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0000000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.00000000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.000000000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0000000000000000001,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0,x0=x0)
    plot_steady_state(sys,symmodel,fact=0.0005,dims=[1,3],x0=x0)
    plot_steady_state(sys,symmodel,dims=[1,3],x0=x0)
    plot_steady_state(sys,symmodel,dims=[1,3],fact=0.0,x0=x0)
end
test_limit_cycle()
end


# function compute_steady_state_Sparse(M)
#     x = [1; (I - M[2:end,2:end]) \ Vector(M[2:end,1])]
#     x =  x/sum(x)
#     a = norm(M * x - x) / norm(x)
#     println(a)
#     return x
# end
#
#
# function compute_steady_state_Sparse_2(M)
#     println(size(M))
#     n = size(M)[1]
#     L = [I-M;ones(1,n)]
#     println(size(L))
#     b = [zeros(n,1);1]
#     x = L\b
#     println(x)
#     return x
# end
#
# function test()
#     # A = sprand(50,50, 1e-2)
#     # M = A./ sum(A, dims=2)
#     # bool = is_Markov_matrix_Sparse(M)
#     # println(bool)
#     #
#     # M = sparse(M')
#     # find an eigenvector x corresponding to eigenvalue λ = 1,
#     # normalized so that the first component is 1
#     A = sprand(1000,1000, 1e-2)
#     A = A ./ sum(A, dims=2)
#     compute_steady_state_Sparse(A)
# end
#
# function test2()
#     I = [1 , 1, 2, 3]
#     J = [2, 3 , 1, 1]
#     V = [0.5 , 0.5, 1.0, 1.0]
#     A = sparse(I,J,V)
#     bool = is_Markov_matrix_Sparse(A)
#     println(bool)
#     compute_steady_state_Sparse_2(A)
#
# end

# test2()
