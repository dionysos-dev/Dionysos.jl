module TestMain

using Test

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PathCompleteFramework

import HybridSystems
import JuMP
import LinearAlgebra as LA
import Clarabel

@testset "PathCompleteFramework" begin
    @testset "edgeList_to_LabDigraph" begin
        edges = [(1, 2, 1.0), (2, 3, 2.0), (1, 3, 3.0)]
        G = PCLF.edgeList_to_LabDigraph(edges)

        @test G isa PCLF.LabDigraph{Float64, Int}
        @test length(G.edges) == 3
        @test G.verts == Set([1, 2, 3])
        @test (1, 2, 1.0) in G.edges
        @test (2, 3, 2.0) in G.edges
        @test (1, 3, 3.0) in G.edges
    end

    @testset "generate_DeBruijn_edges with k = 0" begin
        M = 3
        G = PCLF.generate_DeBruijn_edges(M, 0)

        @test length(G.verts) == 1
        @test 1 in G.verts
        @test length(G.edges) == M
        @test Set(G.edges) == Set([(1, 1, 1), (1, 1, 2), (1, 1, 3)])
    end

    @testset "generate_DeBruijn_edges with k = 1" begin
        M = 2
        G = PCLF.generate_DeBruijn_edges(M, 1)

        expected_verts = Set([(1,), (2,)])
        expected_edges =
            Set([((1,), (1,), 1), ((1,), (2,), 2), ((2,), (1,), 1), ((2,), (2,), 2)])

        @test G.verts == expected_verts
        @test Set(G.edges) == expected_edges
    end

    @testset "generate_DeBruijn_edges with k = 2" begin
        M = 2
        G = PCLF.generate_DeBruijn_edges(M, 2)

        expected_verts = Set([(1, 1), (1, 2), (2, 1), (2, 2)])

        expected_edges = Set([
            ((1, 1), (1, 1), 1),
            ((1, 1), (1, 2), 2),
            ((1, 2), (2, 1), 1),
            ((1, 2), (2, 2), 2),
            ((2, 1), (1, 1), 1),
            ((2, 1), (1, 2), 2),
            ((2, 2), (2, 1), 1),
            ((2, 2), (2, 2), 2),
        ])

        @test G.verts == expected_verts
        @test Set(G.edges) == expected_edges
    end

    @testset "generate_DeBruijn_edges dual orientation" begin
        M = 2
        G = PCLF.generate_DeBruijn_edges(M, 2; dual = true)

        expected_edges = Set([
            ((1, 1), (1, 1), 1),
            ((1, 2), (1, 1), 2),
            ((2, 1), (1, 2), 1),
            ((2, 2), (1, 2), 2),
            ((1, 1), (2, 1), 1),
            ((1, 2), (2, 1), 2),
            ((2, 1), (2, 2), 1),
            ((2, 2), (2, 2), 2),
        ])

        @test Set(G.edges) == expected_edges
    end

    @testset "EllipsoidalPiece sublevel set" begin
        P = Matrix{Float64}(LA.I, 2, 2)
        piece = PCLF.EllipsoidalPiece(P)

        S = PCLF.get_sublevel_set(piece, 1.0)

        @test S !== nothing
    end

    @testset "PolyhedralPiece stores data" begin
        Gmat = [1.0 0.0; 0.0 1.0]
        w = [1.0, 2.0]
        piece = PCLF.PolyhedralPiece(Gmat, w)

        @test piece.G == Gmat
        @test piece.w == w
    end

    @testset "PCLF container" begin
        G = PCLF.generate_DeBruijn_edges(2, 0)
        pieces = Dict{Any, PCLF.AbstractPiece}()
        pieces[1] = PCLF.EllipsoidalPiece(Matrix{Float64}(LA.I, 2, 2))

        pclf = PCLF.PCLF(G, pieces, 0.9)

        @test pclf.graph === G
        @test haskey(pclf.pieces, 1)
        @test pclf.JSRapprox == 0.9
    end

    @testset "compute_quadratic_pieces_pclf on switched system" begin
        A1 = [1.5519 0.4474; 7.6412 7.4716]
        A2 = [0.4750 9.1755; 1.8955 0.1850]
        f = HybridSystems.discreteswitchedsystem([A1, A2])

        G = PCLF.edgeList_to_LabDigraph([
            (1, 3, 2),
            (1, 4, 2),
            (2, 1, 1),
            (2, 2, 1),
            (3, 3, 2),
            (3, 4, 2),
            (4, 1, 1),
            (4, 2, 1),
        ])

        optimizer = JuMP.optimizer_with_attributes(
            Clarabel.Optimizer,
            "max_iter" => 1000,
            "verbose" => false,
        )

        pclf = PCLF.compute_quadratic_pieces_pclf(
            f,
            G,
            optimizer;
            tol = 1e-6,
            maxiter = 100,
            MLF = true,
        )

        @test pclf isa PCLF.PCLF
        @test pclf.JSRapprox >= 0.0
        @test length(pclf.pieces) == 4
        @test Set(keys(pclf.pieces)) == Set([1, 2, 3, 4])

        for v in 1:4
            @test haskey(pclf.pieces, v)
            piece = pclf.pieces[v]
            @test piece isa PCLF.EllipsoidalPiece
            @test size(piece.P) == (2, 2)
            @test isapprox(piece.P, piece.P'; atol = 1e-7)

            evals = LA.eigvals(LA.Symmetric(piece.P))
            @test minimum(real.(evals)) > 0
        end
    end

    @testset "compute_polyhedral_pieces_pclf on switched system" begin
        A1 = [1.5519 0.4474; 7.6412 7.4716]
        A2 = [0.4750 9.1755; 1.8955 0.1850]
        f = HybridSystems.discreteswitchedsystem([A1, A2])

        G = PCLF.generate_DeBruijn_edges(2, 0; dual = false)

        optimizer = JuMP.optimizer_with_attributes(
            Clarabel.Optimizer,
            "max_iter" => 1000,
            "verbose" => false,
        )

        pclf = PCLF.compute_polyhedral_pieces_pclf(
            f,
            G,
            optimizer;
            tol = 1e-6,
            maxiter = 100,
            MLF = true,
            min_w = 1e-5,
        )

        @test pclf isa PCLF.PCLF
        @test pclf.JSRapprox >= 0.0
        @test length(pclf.pieces) == 1
        @test haskey(pclf.pieces, 1)

        piece = pclf.pieces[1]
        @test piece isa PCLF.PolyhedralPiece
        @test size(piece.G) == (2, 2)
        @test length(piece.w) == 2
        @test all(piece.w .>= 1e-5)
    end

    @testset "graph size sanity check" begin
        G = PCLF.generate_DeBruijn_edges(2, 1)
        @test length(G.verts) == 2
        @test length(G.edges) == 4
    end
end

end # module TestMain
