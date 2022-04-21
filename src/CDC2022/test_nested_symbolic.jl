module Test
include("../Dionysos.jl")
using Test, StaticArrays, Plots
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


struct SimpleSystem{N,T,F<:Function}
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

function test()
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = UT.HyperRectangle(SVector(-15.0, -15.0), SVector(-15.0, -15.0))
    hx = [1.0, 1.0]*5.0
    periodic = Int[1,2]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = DO.RectanglularObstacles(X, [obstacle])
    Xdom = DO.GeneralDomainList(hx;elems=d,periodic=periodic,periods=periods,T0=T0,fit=true)

    Udom = DO.CustomList([SVector(0.0,1.0)])
    param = Param(500)
    symmodel = NewNestedSymbolicModel(Xdom, Udom, param)

    # fig = plot(aspect_ratio = 1,legend = false)
    # display(fig)
    s = get_state_by_xpos(symmodel, (2,2), 1)
    s = get_state_by_xpos(symmodel, (2,3), 1)
    println(s)
    plot_automaton(symmodel)

    sys = build_system()
    compute_symbolic_full_domain!(symmodel, sys)
    plot_automaton(symmodel)

    state = get_state_by_coord(symmodel, SVector(4.0,4.0))
    split_node!(symmodel,sys,state)
    plot_automaton(symmodel)

    state = get_state_by_coord(symmodel, SVector(4.0,4.0))
    split_node!(symmodel,sys,state)
    plot_automaton(symmodel)
    update_MC!(symmodel)
    plot_automaton(symmodel)
    plot_steady_state(sys,symmodel)
    plot_shannon_entropy(sys,symmodel)

end

test()



end
