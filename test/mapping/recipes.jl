module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Plots
import LazySets

# The mapping recipes and the pure helpers they lean on. Recipes are asserted through the full
# Plots pipeline (`plot` → `plt.series_list`); see `test/utils/plotting.jl` for why. The helpers
# below it (`merge_rectangles_2d`, `intrect2_to_real_rect`, `project_states_on_dims`) are plain
# functions and are tested directly — they carry the logic that decides *what* gets drawn, so
# pinning them is worth more than pinning a series count.

ENV["GKSwstype"] = "100"

@testset "(grid, pos) recipe" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 2.0))

    plt = plot(grid, (1, 2))
    @test length(plt.series_list) == 1
    s = plt.series_list[1]
    @test s[:seriestype] === :shape
    # Cell (1, 2) is centred on (1, 4) and is one step wide in each direction.
    @test extrema(s[:x]) == (0.5, 1.5)
    @test extrema(s[:y]) == (3.0, 5.0)

    # A cell of a 3-D grid is projected onto `dims` before it is drawn.
    grid3 = MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 4.0))
    s = plot(grid3, (0, 0, 1); dims = [1, 3]).series_list[1]
    @test extrema(s[:x]) == (-0.5, 0.5)
    @test extrema(s[:y]) == (2.0, 6.0)
end

@testset "state-set over grid-mapping recipe" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    MP.cover!(
        m,
        LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(2.0, 2.0)),
        MP.OUTER,
    )
    @test MP.get_n_state(m) == 9

    # Default (`efficient`): the 3×3 block of cells is merged into a single rectangle, which is
    # what keeps a plot of a large abstraction from emitting one series per cell.
    plt = plot(m)
    @test length(plt.series_list) == 1
    @test extrema(plt.series_list[1][:x]) == (-0.5, 2.5)
    @test extrema(plt.series_list[1][:y]) == (-0.5, 2.5)

    # Opting out draws every cell separately.
    @test length(plot((MP.FullStateSet{2}(), m); efficient = false).series_list) == 9

    # Only the first series carries the label, so the legend gets one entry per set.
    @test [s[:label] for s in plot(m; label = "domain").series_list] == ["domain"]
    per_cell = plot((MP.FullStateSet{2}(), m); efficient = false, label = "cells")
    @test [s[:label] for s in per_cell.series_list] == ["cells"; fill("", 8)]

    # A value function forces the per-cell path (each cell gets its own colour) and orders the
    # cells by value, highest first, so cheaper cells are drawn on top.
    vf = q -> Float64(q)
    plt = plot((MP.FullStateSet{2}(), m); value_function = vf)
    @test length(plt.series_list) == 9
    best = MP.get_coord_by_pos(grid, MP.get_pos_by_state(m, 9))
    @test sum(extrema(plt.series_list[1][:x])) / 2 == best[1]
    @test sum(extrema(plt.series_list[1][:y])) / 2 == best[2]

    # An unreachable cell has infinite cost. It must still be drawn (as a gap in the colouring)
    # rather than silently dropped, otherwise the plot would claim the abstraction is smaller
    # than it is.
    @test length(plot((MP.FullStateSet{2}(), m); value_function = q -> Inf).series_list) ==
          9
    part = q -> (q == 1 ? Inf : Float64(q))
    @test length(plot((MP.FullStateSet{2}(), m); value_function = part).series_list) == 9

    # `reducer` decides which value survives when several cells project onto the same pixel.
    @test length(
        plot((MP.FullStateSet{2}(), m); value_function = vf, reducer = max).series_list,
    ) == 9
end

@testset "MappedStateSet recipe" begin
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    MP.cover!(
        m,
        LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(2.0, 2.0)),
        MP.OUTER,
    )

    adjacent = [MP.get_state_by_pos(m, (0, 0)), MP.get_state_by_pos(m, (1, 0))]
    ms = MP.MappedStateSet(MP.stateset_from_states(m, adjacent), m)
    @test length(plot(ms).series_list) == 1

    apart = [MP.get_state_by_pos(m, (0, 0)), MP.get_state_by_pos(m, (2, 0))]
    ms2 = MP.MappedStateSet(MP.stateset_from_states(m, apart), m)
    @test length(plot(ms2).series_list) == 2
end

@testset "merge_rectangles_2d" begin
    @test isempty(MP.merge_rectangles_2d(NTuple{2, Int}[]))

    # A run of consecutive cells in one row collapses to a single rectangle…
    r = MP.merge_rectangles_2d([(0, 0), (1, 0), (2, 0)])
    @test length(r) == 1
    @test r[1].lb == (0, 0) && r[1].ub == (2, 0)

    # …a gap in that run does not.
    @test length(MP.merge_rectangles_2d([(0, 0), (2, 0)])) == 2

    # Rows that span the same columns stack vertically into one rectangle.
    r = MP.merge_rectangles_2d([(0, 0), (1, 0), (0, 1), (1, 1)])
    @test length(r) == 1
    @test r[1].lb == (0, 0) && r[1].ub == (1, 1)

    # Rows that do not agree on their columns cannot stack, and neither can rows that are not
    # vertically adjacent.
    @test length(MP.merge_rectangles_2d([(0, 0), (1, 0), (0, 1)])) == 2
    @test length(MP.merge_rectangles_2d([(0, 0), (1, 0), (0, 2), (1, 2)])) == 2

    # Input order must not matter.
    shuffled = MP.merge_rectangles_2d([(1, 1), (0, 0), (1, 0), (0, 1)])
    @test length(shuffled) == 1
    @test shuffled[1].ub == (1, 1)
end

@testset "intrect2_to_real_rect" begin
    grid = MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(1.0, 2.0, 4.0))
    rect = MP.intrect2_to_real_rect(grid, MP.IntRect2((0, 1), (2, 3)), 1, 3)
    # Bounds run from the outer edge of the first cell to the outer edge of the last.
    @test LazySets.low(rect) == SVector(-0.5, 2.0)
    @test LazySets.high(rect) == SVector(2.5, 14.0)
end

@testset "project_states_on_dims" begin
    # Two layers along the third axis, so several cells project onto the same 2-D pixel.
    grid = MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    m = MP.ExplicitGridMapping(grid)
    MP.cover!(
        m,
        LazySets.Hyperrectangle(;
            low = SVector(0.0, 0.0, 0.0),
            high = SVector(1.0, 1.0, 1.0),
        ),
        MP.OUTER,
    )
    @test MP.get_n_state(m) == 8

    # Without a value function the projection is a plain dedup: one representative cell per
    # pixel, and nothing to colour by.
    proj, posN = MP.project_states_on_dims(m, MP.enum_states(m); dims = [1, 2])
    @test proj === nothing
    @test length(posN) == 4
    @test length(unique(p -> (p[1], p[2]), posN)) == 4

    # With one, the cells sharing a pixel are folded by `reducer`.
    vf = q -> Float64(q)
    proj, posN =
        MP.project_states_on_dims(m, MP.enum_states(m); dims = [1, 2], value_function = vf)
    @test length(posN) == 4
    @test all(v -> v[1] <= 4.0, values(proj))   # `min` keeps the lower of each pair

    proj_max, _ = MP.project_states_on_dims(
        m,
        MP.enum_states(m);
        dims = [1, 2],
        value_function = vf,
        reducer = max,
    )
    @test all(v -> v[1] >= 5.0, values(proj_max))
end

end # module TestMain
