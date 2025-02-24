module TestMain
using StaticArrays, Test, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

sleep(0.1) # used for good printing
println("Started test")

@testset "GeneralDomain" begin
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(10.0, 10.0))
    hx = [2.0, 1.0]
    periodic = [1, 2]
    periods = [10.0, 10.0]

    domain = DO.GeneralDomainList(hx; periodic = periodic, periods = periods)

    DO.add_coord!(domain, SVector(1.2, 3.5))
    @test DO.get_ncells(domain) == 1
    DO.add_coord!(domain, SVector(1.2 + 10.0, 3.5 + 10.0))
    @test DO.get_ncells(domain) == 1
    DO.add_coord!(domain, SVector(15.0, 15.0))
    @test DO.get_ncells(domain) == 2
    DO.add_set!(domain, UT.HyperRectangle(SVector(0.0, 0.0), SVector(15.0, 15.0)), DO.OUTER)
    @test DO.get_ncells(domain) == 50
    DO.remove_pos!(domain, (4, 4))
    @test DO.get_ncells(domain) == 49
    DO.remove_coord!(domain, SVector(11.0, 11.0))
    @test DO.get_ncells(domain) == 48
    pos_iter = DO.enum_pos(domain)
    @test length(pos_iter) == 48
    DO.remove_set!(
        domain,
        UT.HyperRectangle(SVector(4.5, 7.5), SVector(7.5, 12.0)),
        DO.OUTER,
    )
    @test DO.get_ncells(domain) == 36

    fig = plot(; aspect_ratio = :equal)
    plot!(fig, domain)
    @test isa(fig, Plots.Plot{Plots.GRBackend})
end

@testset "GeneralDomain HyperRec representation" begin
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(10.0, 10.0))
    obstacle = UT.HyperRectangle(SVector(4.1, 4.1), SVector(5.1, 5.1))
    d = DO.RectangularObstacles(X, [obstacle])
    hx = [2.0, 1.0]
    periodic = [1, 2]
    periods = [10.0, 10.0]

    domain = DO.GeneralDomainList(
        hx;
        elems = d,
        periodic = periodic,
        periods = periods,
        fit = true,
    )
    @test DO.get_ncells(domain) == 48
end

sleep(0.1) # used for good printing
println("End test")
end # end module
