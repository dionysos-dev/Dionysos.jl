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

function test()
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [3.0, 1.0]*2.0
    periodic = Int[1]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = DO.RectanglularObstacles(X, [obstacle])
    dom = DO.GeneralDomainList(hx;elems=d,periodic=periodic,periods=periods,T0=T0,fit=true)

    Ndomain = NestedDomain(dom)
    println(DO.get_ncells(Ndomain))

    fig = plot(aspect_ratio = 1,legend = false)
    Plots.plot!(Ndomain)
    display(fig)
    cut_pos!(Ndomain, (2,2), 1)
    println(DO.get_ncells(Ndomain))
    cut_pos!(Ndomain, (2,3), 1)
    println(DO.get_ncells(Ndomain))
    cut_pos!(Ndomain, (4,4), 2)
    println(DO.get_ncells(Ndomain))
    fig = plot(aspect_ratio = 1,legend = false)
    Plots.plot!(Ndomain)
    display(fig)
end

test()



end
