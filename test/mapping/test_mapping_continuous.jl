module TestMain
using Test

using StaticArrays

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const MA = DI.Mapping

@testset "MappingContinuousGrid" begin
    orig = SVector(0.0, 0.0)
    h = SVector(0.1, 0.1)
    cont_dom = DO.ContinuousUnboundedDomain(orig)

    x = SVector(0.5, 0.5)

    grid_free = DO.GridFree(orig, h)
    grid_free_dom = DO.DomainList(grid_free)
    mcg_free = MA.MappingContinuousGrid(cont_dom, grid_free_dom)
    @test mcg_free.continuousdomain === cont_dom
    @test mcg_free.griddomain === grid_free_dom
    @test mcg_free.cont2grid(x) === DO.get_rec(grid_free, x)
    xi_free = DO.get_rec(grid_free, x)
    @test mcg_free.grid2cont(xi_free) === UT.get_center(xi_free)

    P = SMatrix{2, 2}(2.0, 1.0, 1.0, 3.0)
    rect = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
    grid_ell = DO.GridEllipsoidalRectangular(orig, h, P)
    grid_ell_dom = DO.DomainList(grid_ell)
    DO.add_set!(grid_ell_dom, rect, DO.INNER)

    mcg_ell = MA.MappingContinuousGrid(cont_dom, grid_ell_dom)
    @test mcg_ell.continuousdomain === cont_dom
    @test mcg_ell.griddomain === grid_ell_dom
    ell_1 = mcg_ell.cont2grid(x)
    ell_2 = DO.get_elem_by_pos(grid_ell, x)
    @test isapprox(ell_1.c, ell_2.c, atol = 1e-6)
    @test isapprox(ell_1.P, ell_2.P, atol = 1e-6)
end

@testset "MappingContinuousEllipsoid" begin
    orig = SVector(0.0, 0.0)
    bound = 3.0
    P = SMatrix{2, 2}(2.0, 1.0, 1.0, 3.0)
    ell = UT.Ellipsoid(P, orig)

    cont_dom = DO.ContinuousUnboundedDomain(orig)
    elli_dom = DO.ContinuousBoundedEllipsoidDomain(orig, bound, Set([ell]))
    mce = MA.MappingContinuousEllipsoid(cont_dom, elli_dom)

    x = SVector(0.5, 0.5)

    @test mce.continuousdomain === cont_dom
    @test mce.ellidomain === elli_dom
    @test mce.cont2elli(x) === nothing
    @test mce.elli2cont(ell) === ell.c
end

end
