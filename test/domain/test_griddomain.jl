module TestMain

using Test
using StaticArrays, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "GridDomain" begin
    orig = SVector(0.0, 0.0)
    h = SVector(1.0, 2.0)
    grid = DO.GridFree(orig, h)
    domain1 = DO.DomainList(grid)

    DO.add_coord!(domain1, SVector(1.2, 3.5))
    @test DO.get_ncells(domain1) == 1
    DO.add_coord!(domain1, SVector(-0.5002, -1.0))
    @test DO.get_ncells(domain1) == 2

    rect = UT.HyperRectangle(SVector(1.0, 0.0), SVector(11.0, 10.0))
    DO.add_set!(domain1, rect, DO.OUTER)
    @test DO.get_ncells(domain1) == 67

    DO.remove_coord!(domain1, SVector(2.0, 2.0))
    @test DO.get_ncells(domain1) == 66
    DO.remove_set!(
        domain1,
        UT.HyperRectangle(SVector(5.0, 5.0), SVector(10000.0, 10000.0)),
        DO.INNER,
    )
    @test DO.get_ncells(domain1) == 48

    pos_iter = DO.enum_pos(domain1)
    @test length(pos_iter) == 48

    domain2 = DO.DomainList(grid)
    union!(domain2, domain1)
    @test DO.get_ncells(domain2) == 48
    DO.remove_set!(
        domain2,
        UT.HyperRectangle(SVector(1.0, 1.0), SVector(2.0, 2.0)),
        DO.OUTER,
    )
    @test DO.get_ncells(domain2) == 45
    DO.add_subset!(
        domain2,
        domain1,
        UT.HyperRectangle(SVector(0.0, 0.0), SVector(5.0, 5.0)),
        DO.INNER,
    )
    @test DO.get_ncells(domain2) == 46

    fig = plot(; aspect_ratio = :equal)
    plot!(fig, domain1)
    @test isa(fig, Plots.Plot{Plots.GRBackend})

    domain3 = DO.DomainList(DO.GridFree(orig, SVector(0.1, 0.1)))
    # First we compile
    @allocated DO.add_set!(domain3, rect, DO.OUTER)
    # Test for the regression detected in https://github.com/dionysos-dev/Dionysos.jl/pull/346#issuecomment-2059069415
    @test (@allocated DO.add_set!(domain3, rect, DO.OUTER)) == 0
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
