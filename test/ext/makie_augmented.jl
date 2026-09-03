module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
using CairoMakie   # loads Makie, activating DionysosMakieExt (headless software renderer)

const PCLFBQ = AB.PCLFBisimulationQuotient

# Build a tiny two-node bisimulation quotient by hand: one square set per node. The augmented
# recipes only read `part_ids`, `states`, and each state's `node`/`set.array`, so this is
# enough to exercise `plot_augmented_bisimulation!` / `plot_augmented_trajectory!` end to end.
function _toy_quotient()
    box1 = LazySets.Hyperrectangle(; low = [0.0, 0.0], high = [1.0, 1.0])
    box2 = LazySets.Hyperrectangle(; low = [1.0, 0.0], high = [2.0, 1.0])
    set1 = LazySets.UnionSetArray([box1])
    set2 = LazySets.UnionSetArray([box2])

    S = typeof(set1)
    U = Tuple{Int}
    slices = Dict{U, Vector{S}}((1,) => [set1], (2,) => [set2])
    bisim = PCLFBQ.PCBisimulationQuotient{S, U}(slices)
    q1 = PCLFBQ.add_state!(bisim, (1,), set1, 0, 1)
    q2 = PCLFBQ.add_state!(bisim, (2,), set2, 0, 1)
    return bisim, q1, q2
end

@testset "Augmented bisimulation Makie recipes" begin
    bisim, q1, q2 = _toy_quotient()
    node_z = Dict((1,) => 1.0, (2,) => 2.0)

    fig = Figure()
    ax = Axis3(fig[1, 1])

    # States: explicit z-map plus each `color_by` branch and the contour outlines.
    @test Dionysos.plot_augmented_bisimulation!(
        ax,
        bisim;
        node_z = node_z,
        color_by = :state,
        show_contours = true,
    ) === ax
    for cb in (:node, :slice, :obs)
        @test Dionysos.plot_augmented_bisimulation!(
            ax,
            bisim;
            node_z = node_z,
            color_by = cb,
        ) === ax
    end
    # Default `node_z` (the extension's internal node-z map) is exercised too.
    @test Dionysos.plot_augmented_bisimulation!(ax, bisim) === ax

    # Batching: two states of different colours become two meshes plus one stroke for both
    # outlines, where drawing them separately gives one mesh and one stroke each. Both must
    # put the same triangles on screen.
    ntris(a) = sum(
        length(Makie.GeometryBasics.faces(p[1][])) for
        p in a.scene.plots if p isa Makie.Mesh
    )
    merged = Axis3(Figure()[1, 1])
    Dionysos.plot_augmented_bisimulation!(
        merged,
        bisim;
        node_z = node_z,
        color_by = :state,
        show_contours = true,
    )
    split = Axis3(Figure()[1, 1])
    Dionysos.plot_augmented_bisimulation!(
        split,
        bisim;
        node_z = node_z,
        color_by = :state,
        show_contours = true,
        merge_plots = false,
    )
    @test length(merged.scene.plots) == 3   # 2 colours + 1 outline stroke
    @test length(split.scene.plots) == 4    # 2 meshes + 2 outline strokes
    @test ntris(merged) == ntris(split)

    # One colour for every state collapses to a single mesh.
    single = Axis3(Figure()[1, 1])
    Dionysos.plot_augmented_bisimulation!(
        single,
        bisim;
        node_z = node_z,
        color_by = :red,
        show_contours = false,
    )
    @test length(single.scene.plots) == 1
    @test ntris(single) == ntris(split)

    # Trajectory augmented onto the node layers.
    X_seq = [[0.5, 0.5], [1.5, 0.5]]
    M_seq = [(0, q1), (0, q2)]
    @test Dionysos.plot_augmented_trajectory!(ax, bisim, X_seq, M_seq; node_z = node_z) ===
          ax

    # A length mismatch of more than one between states and memory is rejected.
    @test_throws ErrorException Dionysos.plot_augmented_trajectory!(
        ax,
        bisim,
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        [(0, q1)];
        node_z = node_z,
    )
end

end # module TestMain
