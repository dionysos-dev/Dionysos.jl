module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Plots

# How a `@recipe` is tested here: run the real Plots pipeline with `plot`, then assert on
# `plt.series_list`. That is what actually exercises the recipe — `plot` expands the recipe
# body, and then keeps expanding whatever each `@series` returns, so a recipe that delegates
# (as most of ours do) is covered together with everything downstream of it. The assertions
# are the recipe's contract: how many series it emits, their coordinates, and which attributes
# it fixes (`:=`) versus merely defaults (`-->`).
#
# Nothing is rendered: `plot` builds a `Plots.Plot`; a backend is only touched on display or
# `savefig`. Headless GR anyway, so this is safe on CI.
ENV["GKSwstype"] = "100"

labels(plt) = [s[:label] for s in plt.series_list]
seriestypes(plt) = [s[:seriestype] for s in plt.series_list]

@testset "DrawPoint recipe" begin
    p = UT.DrawPoint(SVector(1.0, 2.0, 3.0))

    plt = plot(p)
    @test length(plt.series_list) == 1
    s = plt.series_list[1]
    @test s[:seriestype] === :scatter
    @test s[:x] == [1.0]
    @test s[:y] == [2.0]

    # `dims` selects which two coordinates are drawn, as a tuple or a vector.
    for dims in ((1, 3), [1, 3])
        s = plot(p; dims = dims).series_list[1]
        @test s[:x] == [1.0]
        @test s[:y] == [3.0]
    end
    s = plot(p; dims = [3, 2]).series_list[1]
    @test s[:x] == [3.0]
    @test s[:y] == [2.0]

    # The label is a default (`-->`), so the caller wins.
    @test labels(plot(p)) == [""]
    @test labels(plot(p; label = "sample")) == ["sample"]

    # A single index cannot say which plane to draw in.
    @test_throws ArgumentError plot(p; dims = 1)
end

@testset "DrawArrow recipe" begin
    p1, p2 = SVector(0.0, 0.0, 9.0), SVector(1.0, 2.0, 7.0)
    a = UT.DrawArrow(p1, p2)

    plt = plot(a)
    @test length(plt.series_list) == 1
    s = plt.series_list[1]
    @test s[:x] == [0.0, 1.0]
    @test s[:y] == [0.0, 2.0]

    s = plot(a; dims = (1, 3)).series_list[1]
    @test s[:x] == [0.0, 1.0]
    @test s[:y] == [9.0, 7.0]

    # An arrow never carries a label (`label :=`), so it stays out of the legend even when the
    # caller asks — that is what keeps a trajectory's arrows from flooding it.
    @test labels(plot(a; label = "ignored")) == [""]

    # Building from two `DrawPoint`s is the other constructor and must agree.
    a2 = UT.DrawArrow(UT.DrawPoint(p1), UT.DrawPoint(p2))
    @test plot(a2).series_list[1][:x] == [0.0, 1.0]
end

@testset "DrawSegment recipe" begin
    p1, p2 = SVector(0.0, 0.0, 5.0), SVector(1.0, 2.0, 6.0)
    seg = UT.DrawSegment(p1, p2)

    plt = plot(seg)
    @test length(plt.series_list) == 1
    s = plt.series_list[1]
    @test s[:x] == [0.0, 1.0]
    @test s[:y] == [0.0, 2.0]
    @test s[:linestyle] === :dash
    @test labels(plt) == [""]

    # Unlike the other two, this recipe takes `dims` positionally.
    s = plot(seg, [1, 3]).series_list[1]
    @test s[:x] == [0.0, 1.0]
    @test s[:y] == [5.0, 6.0]

    seg2 = UT.DrawSegment(UT.DrawPoint(p1), UT.DrawPoint(p2))
    @test plot(seg2).series_list[1][:y] == [0.0, 2.0]
end

end # module TestMain
