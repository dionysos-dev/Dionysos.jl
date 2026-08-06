module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import LazySets
using Plots

@testset "Set algebra (LazySets-backed)" begin
    A1 = UT.box([-1.0, -1.0], [1.0, 1.0])
    A2 = UT.box([2.0, 2.0], [3.0, 3.0])   # disjoint from A1
    A3 = UT.box([0.5, 0.5], [1.5, 1.5])   # overlaps A1
    B0 = LazySets.Ellipsoid([0.0, 0.0], [1.0 0.0; 0.0 1.0])

    # --- set_union ---
    U = UT.set_union([A1, A2])
    @test U isa LazySets.UnionSetArray
    @test length(U.array) == 2
    @test [0.0, 0.0] ∈ U        # in A1
    @test [2.5, 2.5] ∈ U        # in A2
    @test !([10.0, 10.0] ∈ U)

    # union of a rectangle and an ellipsoid
    Ume = UT.set_union([A2, B0])
    @test [0.0, 0.0] ∈ Ume      # in ellipsoid
    @test [2.5, 2.5] ∈ Ume      # in rectangle
    @test !([1.5, 0.0] ∈ Ume)

    # --- set_minus (A \ B) ---
    M = UT.set_minus(A1, B0)
    @test M isa UT.SetMinus
    @test UT.minus_included(M) === A1
    @test UT.minus_hole(M) === B0
    @test !([0.0, 0.0] ∈ M)     # center is inside the ellipsoid hole
    @test [1.0, 0.5] ∈ M        # inside A1, outside the ellipsoid
    @test !([5.0, 5.0] ∈ M)     # outside A1

    # --- incremental builders ---
    grown = UT.add_set(M, A2)               # enlarge the kept region
    @test [2.5, 2.5] ∈ grown
    @test [1.0, 0.5] ∈ grown
    punched = UT.remove_set(M, A3)          # add A3 to the holes
    @test !([1.2, 1.2] ∈ punched)           # now excluded (in A3)

    # builders are total over bare/empty regions too
    from_empty = UT.add_set(LazySets.EmptySet(2), A1)
    @test [0.0, 0.9] ∈ from_empty
    @test !([5.0, 5.0] ∈ from_empty)

    # --- extractors are total ---
    @test UT.minus_included(A1) === A1
    @test UT.minus_hole(A1) isa LazySets.EmptySet
    e = LazySets.EmptySet(2)
    @test UT.minus_included(e) === e
    @test !([0.0, 0.0] ∈ e)

    # --- dimensions ---
    @test LazySets.dim(U) == 2
    @test LazySets.dim(M) == 2
end

@testset "Inclusion verbs (generic, non-ellipsoid)" begin
    small = UT.box([0.0, 0.0], [1.0, 1.0])
    big = UT.box([-1.0, -1.0], [2.0, 2.0])
    far = UT.box([5.0, 5.0], [6.0, 6.0])

    @test UT.is_included(small, big)
    @test !UT.is_included(big, small)
    @test UT.is_disjoint(small, far)
    @test !UT.is_disjoint(small, big)
end

@testset "Outer bounding box" begin
    H = UT.box([1.0, 2.0], [3.0, 4.0])
    @test UT._outer_box(H) === H                     # a box is its own outer box

    # a set with only a support function goes through box_approximation
    E = LazySets.Ellipsoid([0.0, 0.0], [1.0 0.0; 0.0 1.0])
    BE = UT._outer_box(E)
    @test all(LazySets.low(BE) .≈ [-1.0, -1.0])
    @test all(LazySets.high(BE) .≈ [1.0, 1.0])

    # a union is bounded by the elementwise min/max of its members' boxes
    U = UT.set_union([UT.box([-1.0, -1.0], [1.0, 1.0]), UT.box([2.0, 2.0], [3.0, 3.0])])
    BU = UT._outer_box(U)
    @test all(LazySets.low(BU) .== [-1.0, -1.0])
    @test all(LazySets.high(BU) .== [3.0, 3.0])

    # a set_minus is bounded by its kept region alone
    M = UT.set_minus(H, E)
    @test UT._outer_box(M) === H

    # an empty union has no inferable box
    empty_union = LazySets.UnionSetArray(LazySets.LazySet{Float64}[])
    @test_throws ErrorException UT._outer_box(empty_union)
end

@testset "project_set" begin
    R = UT.box([0.0, 10.0], [2.0, 12.0])
    @test UT.project_set(R, [1, 2]) === R            # dims cover R → identity

    P = UT.project_set(R, [1])                        # drop to the first coordinate
    @test LazySets.dim(P) == 1
    @test [1.0] ∈ P
    @test [3.0] ∉ P
end

@testset "Periodic wrapping" begin
    periodic_dims = UT.SVector(1)
    periods = UT.SVector(1.0)
    start = UT.SVector(0.0)

    R = UT.box([0.8, 0.0], [1.2, 1.0])     # crosses the period boundary
    WR = UT.set_in_period(R, periodic_dims, periods, start)
    @test WR isa LazySets.UnionSetArray
    @test length(WR.array) == 2

    U = UT.set_union([R])
    WU = UT.set_in_period(U, periodic_dims, periods, start)
    @test WU isa LazySets.UnionSetArray
    @test length(WU.array) == 2

    M = UT.set_minus(U, UT.set_union([UT.box([2.0, 2.0], [3.0, 3.0])]))
    WM = UT.set_in_period(M, periodic_dims, periods, start)
    @test WM isa UT.SetMinus
    @test length(UT.minus_included(WM).array) == 2
end

# The recipes are asserted through the full Plots pipeline (`plot` → `plt.series_list`), which
# expands the recipe body and everything it delegates to. See `test/utils/plotting.jl` for why.
@testset "set_union recipe" begin
    U = UT.set_union([UT.box([0.0, 0.0], [1.0, 1.0]), UT.box([3.0, 0.0], [4.0, 1.0])])

    plt = plot(U)
    @test length(plt.series_list) == 2
    # One legend entry for the whole union: only the first part is labelled, so a set drawn in
    # many pieces does not repeat itself in the legend.
    @test [s[:label] for s in plt.series_list] == ["set", ""]
    @test all(s -> s[:seriestype] === :shape, plt.series_list)

    @test [s[:label] for s in plot(U; label = "obstacles").series_list] == ["obstacles", ""]

    # A union living in 3-D is projected by `dims` before it reaches the backend — plotting the
    # unprojected 3-D polytope would need a hull library.
    U3 = UT.set_union([UT.box([0.0, 0.0, 0.0], [1.0, 2.0, 3.0])])
    s = plot(U3; dims = [1, 3]).series_list[1]
    @test maximum(s[:x]) == 1.0
    @test maximum(s[:y]) == 3.0
end

@testset "set_minus recipe" begin
    S = UT.set_minus(UT.box([0.0, 0.0], [4.0, 4.0]), UT.box([1.0, 1.0], [2.0, 2.0]))

    plt = plot(S)
    @test length(plt.series_list) == 2
    # The hole is drawn on top of the enclosing set rather than subtracted from it, and carries
    # no legend entry of its own.
    @test [s[:label] for s in plt.series_list] == ["Set", ""]
    hole = plt.series_list[2]
    @test hole[:seriestype] === :shape
    @test maximum(hole[:x]) == 2.0

    @test [s[:label] for s in plot(S; label = "free space").series_list] == ["free space", ""]
    @test plot(S; hole_alpha = 0.3).series_list[2][:fillalpha] == 0.3

    S3 = UT.set_minus(
        UT.box([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]),
        UT.box([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
    )
    @test length(plot(S3; dims = [1, 3]).series_list) == 2
end

end # module
