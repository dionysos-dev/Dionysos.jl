module Test
using ..Dionysos
using StaticArrays, Plots
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain


function test()
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [3.0, 1.0]*2.0
    periodic = Int[1]
    periods = [30.0,30.0]
    T0 = [0.0,0.0]
    d = DO.RectanglularObstacles(X, [obstacle])
    dom = DO.GeneralDomainList(hx;elems=d,periodic=periodic,periods=periods,T0=T0,fit=true)

    Ndomain = DO.NestedDomain(dom)

    fig = plot(aspect_ratio = 1,legend = false)
    Plots.plot!(Ndomain)
    display(fig)

    DO.cut_pos!(Ndomain, (2,2), 1)

    DO.cut_pos!(Ndomain, (2,3), 1)

    DO.cut_pos!(Ndomain, (4,4), 2)
    fig = plot(aspect_ratio = 1,legend = false)
    DO.plot!(Ndomain)
    display(fig)
end

test()

end
