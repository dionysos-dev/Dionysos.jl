module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

const PCLF = UT.PathCompleteFramework

import HybridSystems
import JuMP
import LinearAlgebra as LA
import Clarabel
import HiGHS

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

    @testset "compute_symmetric_2n_faces_polyhedral_pieces_pclf on switched system" begin
        A1 = [1.5519 0.4474; 7.6412 7.4716]
        A2 = [0.4750 9.1755; 1.8955 0.1850]
        f = HybridSystems.discreteswitchedsystem([A1, A2])

        G = PCLF.generate_DeBruijn_edges(2, 0; dual = false)

        optimizer = JuMP.optimizer_with_attributes(
            Clarabel.Optimizer,
            "max_iter" => 1000,
            "verbose" => false,
        )

        pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
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

        # 2D partition into 4 cones (quadrants)
        e1 = [1.0, 0.0]
        e2 = [0.0, 1.0]

        X1 = hcat(e1, e2)
        X2 = hcat(-e1, e2)
        X3 = hcat(-e1, -e2)
        X4 = hcat(e1, -e2)

        partitions = Dict(1 => [X1, X2, X3, X4])

        pclf = PCLF.compute_polyhedral_pieces_pclf(
            f,
            G,
            optimizer,
            partitions;
            tol = 1e-6,
            maxiter = 100,
            MLF = true,
            min_c = 1e-5,
        )

        @test pclf isa PCLF.PCLF
        @test pclf.JSRapprox >= 0.0
        @test length(pclf.pieces) == 1
        @test haskey(pclf.pieces, 1)

        piece = pclf.pieces[1]
        @test piece isa PCLF.PolyhedralPiece
        @test size(piece.G, 2) == 2
        @test length(piece.w) == size(piece.G, 1)
        @test all(piece.w .== 1.0)
    end
end

# ------------------------------------------------------------------
# Graph shape: completeness and determinism
# ------------------------------------------------------------------
# The two properties are independent and license different uses. Path-completeness represents every
# finite mode word by some path, which is what co-safe synthesis needs; completeness additionally
# gives every node an edge for every mode, which is what verification under arbitrary switching
# needs, since a missing edge removes a move the real adversary has.
@testset "graph completeness and determinism" begin
    for k in 0:2
        G = PCLF.generate_DeBruijn_edges(2, k)
        @test PCLF.is_complete(G, 1:2)
        @test PCLF.is_deterministic(G, 1:2)
    end

    G = PCLF.generate_DeBruijn_edges(2, 1; dual = true)
    @test !PCLF.is_complete(G, 1:2)
    @test !PCLF.is_deterministic(G, 1:2)

    H = PCLF.edgeList_to_LabDigraph([(1, 1, 1), (1, 2, 2), (1, 2, 1), (2, 1, 2)])
    @test !PCLF.is_complete(H, 1:2)          # node 2 has no mode-1 edge
    @test !PCLF.is_deterministic(H, 1:2)     # node 1 branches on mode 1

    for M in (2, 4)
        edges = [(s, mod(s + m - 1, M), m) for s in 0:(M - 1) for m in 1:M]
        R = PCLF.edgeList_to_LabDigraph(edges)
        @test PCLF.is_complete(R, 1:M)
        @test PCLF.is_deterministic(R, 1:M)
    end
end

# ------------------------------------------------------------------
# An absent certificate must be reported, never disguised as a bound
# ------------------------------------------------------------------
# The bisection starts from `b = max_i ||A_i||_inf`, which is *not* known to be feasible. If no trial
# rate is ever feasible, `b` never moves; returning it would report that trivial norm as a computed
# certificate, with the same value for every template -- which is how the failure previously
# masqueraded as a non-monotone bound across template refinements.
@testset "absent certificate reported as Inf, not as the norm bound" begin
    A1 = [1.6 0.0; 0.0 1.6]
    A2 = [0.0 1.6; 1.6 0.0]
    f = HybridSystems.discreteswitchedsystem([A1, A2])
    norm_bound = max(LA.opnorm(A1, Inf), LA.opnorm(A2, Inf))

    G = PCLF.generate_DeBruijn_edges(2, 1)
    nodes = sort(collect(G.verts))
    part = PCLF.conic_partitions_2d(2)
    opt = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 2000,
        "verbose" => false,
    )

    pclf = PCLF.compute_polyhedral_pieces_pclf(
        f,
        G,
        opt,
        Dict(v => part for v in nodes);
        MLF = true,
    )
    @test pclf.JSRapprox == Inf
    @test pclf.JSRapprox != norm_bound
    @test isempty(pclf.pieces)
end

# ------------------------------------------------------------------
# Monotonicity over a nested template family
# ------------------------------------------------------------------
# `conic_partitions_2d` bisects every cone of the previous order, so order p+1 is a strictly richer
# template and the optimal rate cannot increase. Clarabel's default feasibility tolerance (1e-8)
# does not converge at order 3 on this instance; 1e-6 does, and then the sequence is monotone. The
# looser tolerance is part of the fixture, not a workaround for a modelling error.
@testset "bound is non-increasing in nested conic-partition order" begin
    A1 = [0.70 0.10; 0.00 0.65]
    A2 = [0.60 -0.15; 0.10 0.55]
    f = HybridSystems.discreteswitchedsystem([A1, A2])
    G = PCLF.generate_DeBruijn_edges(2, 1)
    nodes = sort(collect(G.verts))
    opt = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 2000,
        "verbose" => false,
        "tol_feas" => 1e-6,
        "tol_gap_abs" => 1e-6,
        "tol_gap_rel" => 1e-6,
    )

    bounds = Float64[]
    for p in 1:3
        part = PCLF.conic_partitions_2d(p)
        pclf = PCLF.compute_polyhedral_pieces_pclf(
            f,
            G,
            opt,
            Dict(v => part for v in nodes);
            MLF = true,
        )
        push!(bounds, pclf.JSRapprox)
    end

    @test all(isfinite, bounds)
    for (prev, next) in zip(bounds, Iterators.drop(bounds, 1))
        @test next <= prev + 1e-3
    end
end

@testset "conic partitions refine" begin
    for p in 1:3
        @test length(PCLF.conic_partitions_2d(p)) == 4 * 2^(p - 1)
    end
end

@testset "mode_matrices accepts both resetmap shapes" begin
    A1, A2 = [0.5 0.0; 0.0 0.4], [0.3 0.1; 0.0 0.2]
    f = HybridSystems.discreteswitchedsystem([A1, A2])

    # `discreteswitchedsystem` wraps a bare matrix in a `LinearMap`, so this exercises the
    # `:A in fieldnames` branch; the bare-matrix branch is exercised by passing them through.
    @test UT.mode_matrices(f) == [A1, A2]
    @test UT.mode_matrices(f) isa Vector{Matrix{Float64}}
end

@testset "check_pclf reproduces a known contraction rate" begin
    # x ↦ 0.5x under V(x) = ‖x‖_∞ contracts at exactly 0.5.
    f = HybridSystems.discreteswitchedsystem([Matrix(0.5 * LA.I, 2, 2)])
    graph = PCLF.edgeList_to_LabDigraph([(1, 1, 1)])
    pieces = Dict{Any, PCLF.AbstractPiece}(
        1 => PCLF.PolyhedralPiece(Matrix(1.0LA.I, 2, 2), ones(2)),
    )

    ok = PCLF.check_pclf(PCLF.PCLF(graph, pieces, 0.5), f)
    @test ok.rate ≈ 0.5 atol = 1e-6
    @test isempty(ok.violated)

    # A rate the system does not attain must be refuted.
    @test length(PCLF.check_pclf(PCLF.PCLF(graph, pieces, 0.2), f).violated) == 1
end

@testset "check_pclf normalises the homogeneity degree" begin
    # The same system certified quadratically: the value ratio is 0.25 but the gauge rate is 0.5.
    # Without the degree correction a quadratic certificate would read as ρ² and look better than
    # it is, and would not be comparable with a polyhedral one on the same system.
    f = HybridSystems.discreteswitchedsystem([Matrix(0.5 * LA.I, 2, 2)])
    graph = PCLF.edgeList_to_LabDigraph([(1, 1, 1)])
    quad = Dict{Any, PCLF.AbstractPiece}(1 => PCLF.EllipsoidalPiece(Matrix(1.0LA.I, 2, 2)))

    @test PCLF.piece_degree(PCLF.EllipsoidalPiece(Matrix(1.0LA.I, 2, 2))) == 2
    @test PCLF.piece_degree(PCLF.PolyhedralPiece(Matrix(1.0LA.I, 2, 2), ones(2))) == 1
    @test PCLF.check_pclf(PCLF.PCLF(graph, quad, 0.5), f).rate ≈ 0.5 atol = 1e-6
end

@testset "certify_pclf is exact where check_pclf only samples" begin
    # Example 3.1 of Gol, Ding, Lazar & Belta with their published certificate. The exact rate is
    # 0.940008, so their printed 0.94 is that supremum rounded: a sample lands just under it and
    # the linear program lands just over. Both must agree to six decimals on the value itself.
    A1 = [-0.65 0.32; -0.42 -0.92]
    A2 = [0.65 0.32; -0.42 -0.92]
    f = HybridSystems.discreteswitchedsystem([A1, A2])
    L = [-0.0625 1.0; 0.6815 1.0; 0.9947 0.6868; 0.9947 -0.0678]

    graph = PCLF.generate_DeBruijn_edges(2, 0)
    pclf = PCLF.PCLF(
        graph,
        Dict{Any, PCLF.AbstractPiece}(1 => PCLF.PolyhedralPiece(L, ones(size(L, 1)))),
        0.94,
    )

    exact = PCLF.certify_pclf(pclf, f, HiGHS.Optimizer)
    sampled = PCLF.check_pclf(pclf, f)

    @test exact.rate ≈ 0.940008 atol = 1e-5
    # A sample of a supremum can only fall short of it.
    @test sampled.rate <= exact.rate + 1e-9

    # Refutation must be unanimous when the claimed rate is far from the true one.
    bad = PCLF.PCLF(graph, pclf.pieces, 0.5)
    @test length(PCLF.certify_pclf(bad, f, HiGHS.Optimizer).violated) == 2
    @test length(PCLF.check_pclf(bad, f).violated) == 2
end

@testset "checkers reject what they cannot decide" begin
    f = HybridSystems.discreteswitchedsystem([Matrix(0.5 * LA.I, 2, 2)])
    graph = PCLF.edgeList_to_LabDigraph([(1, 1, 1)])

    # The failure the solvers now signal with `JSRapprox = Inf` carries no pieces at all.
    empty_pclf = PCLF.PCLF(graph, Dict{Any, PCLF.AbstractPiece}(), Inf)
    @test_throws ErrorException PCLF.check_pclf(empty_pclf, f)

    # `certify_pclf` is exact only for polyhedral pieces and must say so rather than approximate.
    quad = Dict{Any, PCLF.AbstractPiece}(1 => PCLF.EllipsoidalPiece(Matrix(1.0LA.I, 2, 2)))
    @test_throws ErrorException PCLF.certify_pclf(
        PCLF.PCLF(graph, quad, 0.5),
        f,
        HiGHS.Optimizer,
    )
end

end # module TestMain
