module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Plots
import LazySets

const PCQ = AB.PCLFBisimulationQuotient

# The quotient recipes drive the figures of the PCLF case study. The solver that normally
# produces a quotient needs an SDP solver and a few seconds, so the quotient here is assembled by
# hand: the recipes only read `states`/`slices`, and building those directly keeps this suite
# fast. Assertions run through the full Plots pipeline; see `test/utils/plotting.jl` for why.

ENV["GKSwstype"] = "100"

sls(boxes...) = UT.semilinear_set([
    LazySets.Hyperrectangle(; low = lo, high = hi) for (lo, hi) in boxes
])

# Two nodes, two slices each. State 3 is deliberately made of two pieces, which is the case that
# distinguishes "one legend entry per abstract state" from "one per polytope".
const S11 = sls(([0.0, 0.0], [1.0, 1.0]))
const S12 = sls(([1.0, 0.0], [2.0, 1.0]))
const S21 = sls(([0.0, 1.0], [1.0, 2.0]), ([2.0, 2.0], [3.0, 3.0]))

function build_quotient()
    S = typeof(S11)
    slices = Dict{Int, Vector{S}}(1 => [S11, S12], 2 => [S21])
    T = PCQ.PCBisimulationQuotient{S, Int}(slices)
    PCQ.add_state!(T, 1, S11, 1, 1)     # node 1, obs 1, slice 1
    PCQ.add_state!(T, 1, S12, 0, 2)     # node 1, obs 0, slice 2
    PCQ.add_state!(T, 2, S21, 1, 1)     # node 2, obs 1, slice 1 — two polytopes
    return T
end

@testset "PCAbstractState recipe" begin
    T = build_quotient()
    q1, q3 = T.states[1], T.states[3]

    plt = plot(q1)
    @test length(plt.series_list) == 1
    # The label names the state and where it came from, so a figure of many states stays readable.
    @test occursin("q1", plt.series_list[1][:label])
    @test occursin("node=1", plt.series_list[1][:label])
    @test plt.series_list[1][:fillcolor] == Plots.plot_color(:blue)

    @test plot(q1; color = :red).series_list[1][:fillcolor] == Plots.plot_color(:red)
    @test [s[:label] for s in plot(q1; show_label = false).series_list] == [""]

    # A state made of several polytopes is still one legend entry.
    plt = plot(q3)
    @test length(plt.series_list) == 2
    @test count(s -> !isempty(s[:label]), plt.series_list) == 1
end

@testset "PCBisimulationQuotient recipe: states" begin
    T = build_quotient()

    # One series per polytope, the unmerged layout. `merge_series = false` is what a caller
    # passes when the drawing order of individual cells matters.
    @test length(plot(T; merge_series = false).series_list) == 4

    # Filters select a subset. Each is an independent predicate.
    unmerged(; kw...) = plot(T; merge_series = false, kw...).series_list
    @test length(unmerged(; node = 1)) == 2
    @test length(unmerged(; node = 2)) == 2   # one state, two polytopes
    @test length(unmerged(; slice = 1)) == 3
    @test length(unmerged(; obs = 0)) == 1
    @test length(unmerged(; state_ids = [1])) == 1
    @test length(unmerged(; state_ids = [1, 2])) == 2
    @test isempty(unmerged(; node = 1, obs = 0, slice = 1))

    # Merged is the default: cells sharing a colour become one NaN-separated series, so the
    # count follows the number of distinct colours rather than the number of polytopes. Three
    # states coloured by state give three; the two obs-1 states share one when coloured by obs.
    @test length(plot(T).series_list) == 3
    @test length(plot(T; by = :obs).series_list) == 2
    @test length(plot(T; user_color = :black).series_list) == 1
    @test length(plot(T; node = 2).series_list) == 1   # two polytopes, one colour
    @test isempty(plot(T; node = 1, obs = 0, slice = 1).series_list)

    # Merging must not drop geometry: the same polytopes are still there, as NaN-separated runs.
    @test count(isnan, plot(T; node = 2).series_list[1][:x]) == 2

    # Labels are off by default — with one entry per abstract state the legend would swamp the
    # figure. Switching them on needs a series per state, so it turns merging off by itself and
    # still gives one entry per state, not one per polytope.
    @test all(s -> isempty(s[:label]), plot(T).series_list)
    @test count(s -> !isempty(s[:label]), plot(T; show_labels = true).series_list) == 3
    @test length(plot(T; show_labels = true).series_list) == 4

    # `by` chooses what the colouring means. Colouring by observation must give the two obs-1
    # states the same colour, which colouring by state must not.
    colours(; kw...) = [s[:fillcolor] for s in unmerged(; kw...)]
    @test colours(; by = :obs)[1] == colours(; by = :obs)[3]
    @test colours(; by = :state)[1] != colours(; by = :state)[2]
    @test length(unique(colours(; by = :slice))) == 2
    @test length(unique(colours(; by = :node))) <= 2
    # An explicit colour overrides every scheme.
    @test length(unique(colours(; user_color = :black))) == 1

    # Cell outlines are hidden unless asked for, so adjacent cells read as one region.
    @test all(s -> s[:linealpha] == 0.0, plot(T).series_list)
    @test all(s -> s[:linealpha] == 1.0, plot(T; show_contours = true).series_list)

    @test all(s -> s[:fillalpha] == 0.5, plot(T; fillalpha = 0.5).series_list)
end

@testset "PCBisimulationQuotient recipe: slices" begin
    T = build_quotient()

    # The slices are the partition the quotient was refined over, independent of its states.
    slices(; kw...) = plot(T; what = :slices, merge_series = false, kw...).series_list
    @test length(slices()) == 4
    @test length(slices(; node = 1)) == 2
    @test length(slices(; node = 2)) == 2
    @test length(slices(; slice = 2)) == 1

    # Merged, slices sharing a colour (same index, different node) become one series.
    @test length(plot(T; what = :slices).series_list) == 2
    @test length(plot(T; what = :slices, node = 1).series_list) == 2
    @test length(plot(T; what = :slices, node = 2).series_list) == 1

    @test all(s -> isempty(s[:label]), plot(T; what = :slices).series_list)
    labelled = plot(T; what = :slices, show_labels = true)
    @test count(s -> !isempty(s[:label]), labelled.series_list) == 3
    @test any(s -> occursin("slice=", s[:label]), labelled.series_list)

    # Neither view is drawn when neither is asked for.
    @test isempty(plot(T; what = :none).series_list)
end

end # module TestMain
