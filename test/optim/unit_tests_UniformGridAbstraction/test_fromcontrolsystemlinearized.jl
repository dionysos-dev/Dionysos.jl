module TestMain

using Test
using StaticArrays, MathematicalSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "FromControlSystem" begin
    lb = SVector(0.0, 0.0)
    ub = SVector(10.0, 11.0)
    x0 = SVector(0.0, 0.0)
    h = SVector(1.0, 2.0)
    Xgrid = DO.GridFree(x0, h)
    Xfull = DO.DomainList(Xgrid)
    DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)

    lb = SVector(-1.0)
    ub = SVector(1.0)
    u0 = SVector(0.0)
    h = SVector(0.5)
    Ugrid = DO.GridFree(u0, h)
    Ufull = DO.DomainList(Ugrid)
    DO.add_set!(Ufull, UT.HyperRectangle(lb, ub), DO.OUTER)

    tstep = 0.5
    # F_sys(x, u) = (1.0-cos(x(2)), -x(1) + u(1)
    # L_growthbound(u) = (0.0 1.0; 1.0 0.0)
    F_sys(x, u) = SVector(u[1], -cos(x[1]))
    DF_sys(x, u) = SMatrix{2, 2}(0.0, sin(x[1]), 0.0, 0.0)
    bound_DF(u) = 1.0
    # DDF_1 = [0.0 0.0; 0.0 0.0]
    # DDF_2 = [cos(x[1]) 0.0; 0.0 0.0]
    bound_DDF(u) = 1.0

    concrete_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_sys,
        2,
        1,
        nothing,
        nothing,
    )
    continuous_approx =
        ST.ContinuousTimeLinearized(concrete_system, DF_sys, bound_DF, bound_DDF)
    discrete_approx = ST.discretize(continuous_approx, tstep)

    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
    SY.compute_abstract_system_from_concrete_system!(symmodel, discrete_approx)

    @test SY.ntransitions(symmodel.autom) == 2175

    xpos = (1, 2)
    x = DO.get_coord_by_pos(Xgrid, xpos)
    source = SY.get_state_by_xpos(symmodel, xpos)
    symbol = 1
    u = SY.get_concrete_input(symmodel, symbol)

    Xsimple = DO.DomainList(Xgrid)
    DO.add_pos!(Xsimple, xpos)
    Ysimple = DO.DomainList(Xgrid)
    targetlist = Int[]
    SY.compute_post!(targetlist, symmodel.autom, source, symbol)
    for target in targetlist
        DO.add_pos!(Ysimple, SY.get_xpos_by_state(symmodel, target))
    end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
