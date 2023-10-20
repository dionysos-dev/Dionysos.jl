module TestMain

using Test
using StaticArrays
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

    tstep = 5.0
    nsys = 3
    ngrowthbound = 3
    # F_sys(x, u) = [1.0-cos(x[2]), -x[1] + u[1]]
    # L_growthbound(u) = [0.0 1.0; 1.0 0.0]
    F_sys(x, u) = SVector(u[1], -cos(x[1]))
    L_growthbound(u) = SMatrix{2, 2}(0.0, 1.0, 0.0, 0.0)
    sysnoise = SVector(1.0, 1.0) * 0.1
    measnoise = SVector(1.0, 1.0) * 0.0

    contsys = ST.NewControlSystemGrowthRK4(
        tstep,
        F_sys,
        L_growthbound,
        sysnoise,
        measnoise,
        nsys,
        ngrowthbound,
    )
    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)
    SY.compute_symmodel_from_controlsystem!(symmodel, contsys)
    @test SY.ntransitions(symmodel.autom) == 1145

    xpos = (1, 2)
    upos = (1,)
    x = DO.get_coord_by_pos(Xgrid, xpos)
    u = DO.get_coord_by_pos(Ugrid, upos)
    source = SY.get_state_by_xpos(symmodel, xpos)
    symbol = SY.get_symbol_by_upos(symmodel, upos)

    Xsimple = DO.DomainList(Xgrid)
    DO.add_pos!(Xsimple, xpos)
    Usimple = DO.DomainList(Ugrid)
    DO.add_pos!(Usimple, upos)
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
