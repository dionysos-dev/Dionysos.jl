module Test
using Test, StaticArrays, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

sleep(0.1) # used for good printing
println("Started test")

@testset "NestedDomain" begin
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))
    hx = [3.0, 1.0] * 2.0
    periodic = Int[1]
    periods = [30.0, 30.0]
    T0 = [0.0, 0.0]
    d = DO.RectangularObstacles(X, [obstacle])
    dom = DO.GeneralDomainList(
        hx;
        elems = d,
        periodic = periodic,
        periods = periods,
        T0 = T0,
        fit = true,
    )

    Ndomain = DO.NestedDomain(dom)

    x = SVector(12.0, 4.5)
    @test DO.get_levels(Ndomain) == 1
    @test DO.get_ncells(Ndomain) == 67
    @test DO.get_depth(Ndomain, x) == 1

    DO.cut_pos!(Ndomain, (2, 2), 1)
    @test DO.get_levels(Ndomain) == 2
    @test DO.get_ncells(Ndomain) == 70
    @test DO.get_depth(Ndomain, x) == 2

    DO.cut_pos!(Ndomain, (2, 3), 1)
    @test DO.get_levels(Ndomain) == 2
    @test DO.get_ncells(Ndomain) == 73
    @test DO.get_depth(Ndomain, x) == 2

    DO.cut_pos!(Ndomain, (4, 4), 2)
    @test DO.get_levels(Ndomain) == 3
    @test DO.get_ncells(Ndomain) == 76
    @test DO.get_depth(Ndomain, x) == 3

    fig = plot(; aspect_ratio = :equal, legend = false)
    plot!(fig, Ndomain)
    @test isa(fig, Plots.Plot{Plots.GRBackend})
end

end
