module TestMain

using Test
using StaticArrays
import LinearAlgebra as LA
import LazySets
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping

println("Started grid tests")

# Helper: check point is inside a box (inclusive)
_in_rec(rect, x) = x ∈ rect

@testset "Grid API" begin
    @testset "GridFree basics" begin
        h = SVector(1.0, 2.0)
        g = MP.GridFree(h; zero_origin = true)

        @test MP.get_origin(g) == SVector(0.0, 0.0)
        @test MP.get_h(g) == h
        @test MP.get_dim(g) == 2

        pos = (3, -2)
        x = MP.get_coord_by_pos(g, pos)
        @test x == SVector(3.0, -4.0)
        @test MP.get_pos_by_coord(g, x) == pos

        @test MP.get_pos_by_coord(g, SVector(0.49, 0.99*2.0)) == (0, 1)
        @test MP.get_pos_by_coord(g, SVector(0.51, 0.01)) == (1, 0)

        r = MP.get_rec(g, (0, 0))
        @test LazySets.low(r) == SVector(-0.5, -1.0)
        @test LazySets.high(r) == SVector(0.5, 1.0)

        @test MP.get_elem_by_pos(g, (0, 0)) == r
        @test MP.get_elem_by_coord(g, SVector(0.1, -0.1)) == MP.get_rec(g, (0, 0))

        @test isapprox(MP.get_volume(g), 2.0; atol = 1e-12)

        @test MP.is_state_cover(g) == false
    end

    @testset "Position limits" begin
        h = SVector(1.0, 1.0)
        g = MP.GridFree(h; zero_origin = true)

        rect = UT.box(SVector(-0.51, -0.51), SVector(0.51, 0.51))

        inner = MP.get_pos_lims_inner(g, rect)
        @test first.(inner) == (0, 0)
        @test last.(inner) == (0, 0)

        outer = MP.get_pos_lims_outer(g, rect)
        @test all(i -> first(outer[i]) <= first(inner[i]), 1:2)
        @test all(i -> last(outer[i]) >= last(inner[i]), 1:2)

        @test MP.get_pos_lims(g, rect, MP.INNER) == inner
        @test MP.get_pos_lims(g, rect, MP.OUTER) == outer
    end

    @testset "get_pos_from_set box" begin
        h = SVector(1.0, 1.0)
        g = MP.GridFree(h; zero_origin = true)
        rect = UT.box(SVector(-1.1, -0.1), SVector(1.1, 0.1))
        poss = collect(MP.get_pos_from_set(g, rect, MP.OUTER))
        @test !isempty(poss)
        @test (0, 0) in poss
    end

    @testset "get_pos_from_set set_union" begin
        h = SVector(1.0, 1.0)
        g = MP.GridFree(h; zero_origin = true)
        r1 = UT.box(SVector(-0.1, -0.1), SVector(0.1, 0.1))
        r2 = UT.box(SVector(2.9, -0.1), SVector(3.1, 0.1))
        U = UT.set_union([r1, r2])
        poss = collect(MP.get_pos_from_set(g, U, MP.OUTER))
        @test (0, 0) in poss
        @test (3, 0) in poss
    end

    @testset "get_pos_from_set set_minus" begin
        h = SVector(1.0, 1.0)
        g = MP.GridFree(h; zero_origin = true)
        A = UT.box(SVector(-0.6, -0.6), SVector(0.6, 0.6))
        B = UT.box(SVector(-0.51, -0.51), SVector(0.51, 0.51))
        S = UT.set_minus(A, B)
        poss = collect(MP.get_pos_from_set(g, S, MP.OUTER))
        @test (0, 0) ∉ poss
        @test !isempty(poss)
    end

    @testset "set_minus respects inverse inclusion mode" begin
        h = SVector(1.0, 1.0)
        g = MP.GridFree(h; zero_origin = true)
        A = UT.box(SVector(-0.6, -0.6), SVector(0.6, 0.6))
        Btiny = UT.box(SVector(-0.1, -0.1), SVector(0.1, 0.1))
        S1 = UT.set_minus(A, Btiny)
        poss1 = collect(MP.get_pos_from_set(g, S1, MP.OUTER))
        @test (0, 0) in poss1
    end

    @testset "GridEllipsoidalRectangular" begin
        orig = SVector(0.0, 0.0)
        h = SVector(1.0, 1.0)
        P = @SMatrix [1.0 0.0; 0.0 1.0]

        g = MP.GridEllipsoidalRectangular(orig, h, P)

        @test MP.get_origin(g) == orig
        @test MP.get_h(g) == h
        @test MP.is_state_cover(g)

        e = MP.get_elem_by_pos(g, (2, -1))
        @test e isa UT.Ellipsoid

        x0 = MP.get_coord_by_pos(g, (0, 0))
        poss = MP.get_all_pos_by_coord(g, x0)

        @test (0, 0) in poss

        x = SVector(0.6, 0.0)
        poss2 = MP.get_all_pos_by_coord(g, x)

        @test !isempty(poss2)
    end
end

end
