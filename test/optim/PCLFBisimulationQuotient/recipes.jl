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

    # Three states, four polytopes between them: one series per polytope.
    @test length(plot(T).series_list) == 4

    # Filters select a subset. Each is an independent predicate.
    @test length(plot(T; node = 1).series_list) == 2
    @test length(plot(T; node = 2).series_list) == 2   # one state, two polytopes
    @test length(plot(T; slice = 1).series_list) == 3
    @test length(plot(T; obs = 0).series_list) == 1
    @test length(plot(T; state_ids = [1]).series_list) == 1
    @test length(plot(T; state_ids = [1, 2]).series_list) == 2
    @test isempty(plot(T; node = 1, obs = 0, slice = 1).series_list)

    # Labels are off by default — with one entry per abstract state the legend would swamp the
    # figure — and switching them on gives one per state, not one per polytope.
    @test all(s -> isempty(s[:label]), plot(T).series_list)
    @test count(s -> !isempty(s[:label]), plot(T; show_labels = true).series_list) == 3

    # `by` chooses what the colouring means. Colouring by observation must give the two obs-1
    # states the same colour, which colouring by state must not.
    colours(; kw...) = [s[:fillcolor] for s in plot(T; kw...).series_list]
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
    @test length(plot(T; what = :slices).series_list) == 4
    @test length(plot(T; what = :slices, node = 1).series_list) == 2
    @test length(plot(T; what = :slices, node = 2).series_list) == 2
    @test length(plot(T; what = :slices, slice = 2).series_list) == 1

    @test all(s -> isempty(s[:label]), plot(T; what = :slices).series_list)
    labelled = plot(T; what = :slices, show_labels = true)
    @test count(s -> !isempty(s[:label]), labelled.series_list) == 3
    @test any(s -> occursin("slice=", s[:label]), labelled.series_list)

    # Neither view is drawn when neither is asked for.
    @test isempty(plot(T; what = :none).series_list)
end

end # module TestMain
