module Test
include("../Dionysos.jl")
using Dionysos, StaticArrays, Plots
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic


function sys_map(x::SVector{N,T}, u, tstep) where {N,T}
    return x+tstep*u
end

function build_system()
    tstep = 0.8
    measnoise = SVector(0.0, 0.0)
    return ST.SimpleSystem(tstep, measnoise, sys_map, nothing)
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
    param = SY.Param(500)
    symmodel = SY.NewNestedSymbolicModel(Xdom, Udom, param)

    s = SY.get_state_by_xpos(symmodel, (2,2), 1)
    s = SY.get_state_by_xpos(symmodel, (2,3), 1)
    SY.plot_automaton(symmodel)

    sys = build_system()
    SY.compute_symbolic_full_domain!(symmodel, sys)
    SY.plot_automaton(symmodel)

    state = SY.get_state_by_coord(symmodel, SVector(4.0,4.0))
    SY.split_node!(symmodel,sys,state)
    SY.plot_automaton(symmodel)

    state = SY.get_state_by_coord(symmodel, SVector(4.0,4.0))
    SY.split_node!(symmodel,sys,state)
    SY.plot_automaton(symmodel)
    SY.update_MC!(symmodel)
    SY.plot_automaton(symmodel)
    SY.plot_steady_state(sys,symmodel)
    SY.plot_shannon_entropy(sys,symmodel)
end

test()



end
