module TestCarving

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import Random

const GRID = MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.1))
const BOX = LazySets.Hyperrectangle(; low = SVector(-0.5, -0.5), high = SVector(0.5, 0.5))

@testset "CellUnion basics" begin
    S = MP.CellUnion(GRID, [(0, 0), (1, 0), (1, 0)])
    @test length(S) == 2 # deduplicated
    @test !isempty(S)
    @test (1, 0) in S
    @test (2, 0) ∉ S
    @test Set(S) == Set([(0, 0), (1, 0)]) # iterable of positions

    # Coordinate membership quantizes to the containing cell.
    @test SVector(0.08, 0.01) in S
    @test SVector(0.08, 0.06) ∉ S # cell (1, 1)

    @test LazySets.dim(S) == 2
    box = UT._outer_box(S)
    @test LazySets.low(box) ≈ [-0.05, -0.05]
    @test LazySets.high(box) ≈ [0.15, 0.05]
end

@testset "exact discretization under every inclusion mode" begin
    S = MP.CellUnion(GRID, [(0, 0), (1, 0)])
    for mode in (MP.INNER, MP.OUTER, MP.CENTER)
        # Exactly the member cells — in particular no face-adjacent neighbor
        # under OUTER, which is what makes cell-aligned holes and targets
        # exact without any ε resizing.
        @test Set(MP.get_pos_from_set(GRID, S, mode)) == Set([(0, 0), (1, 0)])
    end

    # A CellUnion built on a different grid must not be silently re-quantized.
    other = MP.GridFree(SVector(0.05, 0.0), SVector(0.1, 0.1))
    @test_throws ErrorException MP.get_pos_from_set(other, S, MP.OUTER)
end

@testset "cell-aligned hole in a set_minus domain" begin
    hole = MP.CellUnion(GRID, [(0, 0), (1, 0)])
    X = UT.set_minus(BOX, hole)
    @test X isa UT.SetMinus

    kept = Set(MP.get_pos_from_set(GRID, X, MP.INNER))
    @test (0, 0) ∉ kept && (1, 0) ∉ kept
    # The face-adjacent neighbors survive — the former ε-shrink regression.
    for neighbor in ((-1, 0), (2, 0), (0, 1), (0, -1), (1, 1), (1, -1))
        @test neighbor in kept
    end

    # Concrete membership is consistent with the discretization.
    @test SVector(0.3, 0.3) in X
    @test SVector(0.02, -0.01) ∉ X # inside the hole
end

@testset "cells_where and mapping-level recovery" begin
    S = MP.cells_where(GRID, BOX) do pos
        x = MP.get_coord_by_pos(GRID, pos)
        return x[1] > 0.24 && x[2] > 0.24
    end
    @test Set(S) == Set((i, j) for i in 3:4, j in 3:4)

    m = MP.ExplicitGridMapping{2, Float64}(GRID)
    MP.cover!(m, BOX, MP.INNER)
    states = MP.get_states_from_set(m, S, MP.INNER)
    @test Set(MP.get_pos_by_state(m, q) for q in states) == Set(S)
end

@testset "image_blocked_cells is sound (adversarial sampling)" begin
    # Nonlinear image map with per-axis gradient bounds 2 and 1.
    g(x) = SVector(2 * x[1] + 0.1 * sin(10 * x[2]), x[2])
    grad = SVector(2.0 + 1.0, 1.0) # |∂g₁/∂x₂| ≤ 1 from the sine term
    obstacle = LazySets.Hyperrectangle(; low = SVector(0.1, 0.1), high = SVector(0.3, 0.3))

    blocked = MP.image_blocked_cells(g, grad, GRID, BOX, obstacle)
    @test !isempty(blocked)

    rng = Random.MersenneTwister(1)
    n_checked = 0
    n_violations = 0
    for _ in 1:20000
        x = SVector(0.5, 0.5) .* (2 .* rand(rng, 2) .- 1)
        pos = MP.get_pos_by_coord(GRID, x)
        pos ∉ blocked || continue
        n_checked += 1
        # A kept cell never maps into the obstacle: that is the certificate.
        g(x) ∈ obstacle && (n_violations += 1)
    end
    @test n_checked > 0
    @test n_violations == 0
end

@testset "swept_input_filter blocks crossing segments" begin
    blocked = MP.CellUnion(GRID, [(1, 0)])
    filter = MP.swept_input_filter(GRID, 1.0, blocked)

    # A two-cell step across the blocked cell is rejected; going around it or
    # stopping short is not.
    @test !filter(SVector(0.0, 0.0), SVector(0.2, 0.0))
    @test !filter(SVector(0.0, 0.0), SVector(0.1, 0.0))
    @test filter(SVector(0.0, 0.0), SVector(-0.1, 0.0))
    @test filter(SVector(0.0, 0.1), SVector(0.2, 0.0))
end

end # module
