module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

const PCLF = UT.PathCompleteFramework

import HybridSystems
import LinearAlgebra as LA

# The De Bruijn graphs these tests read, labelled the way Dionysos labels them: a node is the tuple
# of the last `k` modes, oldest first. Dionysos no longer builds certificates, so it no longer
# builds their graphs either; PathCompleteCertificates.jl does both (see test/ext/).
function de_bruijn(M, k; dual = false)
    k == 0 && return PCLF.edgeList_to_LabDigraph([(1, 1, m) for m in 1:M])
    nodes = vec(collect(Iterators.product(ntuple(_ -> 1:M, k)...)))
    edges = [(s, (s[2:end]..., m), m) for s in nodes for m in 1:M]
    dual && (edges = [(d, s, m) for (s, d, m) in edges])
    return PCLF.edgeList_to_LabDigraph(edges)
end

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
        @test PCLF.piece_value(piece, [0.5, -3.0]) == 1.5
    end

    @testset "PCLF container" begin
        G = de_bruijn(2, 0)
        pieces = Dict{Any, PCLF.AbstractPiece}()
        pieces[1] = PCLF.EllipsoidalPiece(Matrix{Float64}(LA.I, 2, 2))

        pclf = PCLF.PCLF(G, pieces, 0.9)

        @test pclf.graph === G
        @test haskey(pclf.pieces, 1)
        @test pclf.JSRapprox == 0.9
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
        G = de_bruijn(2, k)
        @test PCLF.is_complete(G, 1:2)
        @test PCLF.is_deterministic(G, 1:2)
    end

    G = de_bruijn(2, 1; dual = true)
    @test !PCLF.is_complete(G, 1:2)
    @test !PCLF.is_deterministic(G, 1:2)

    H = PCLF.edgeList_to_LabDigraph([(1, 1, 1), (1, 2, 2), (1, 2, 1), (2, 1, 2)])
    @test !PCLF.is_complete(H, 1:2)          # node 2 has no mode-1 edge
    @test !PCLF.is_deterministic(H, 1:2)     # node 1 branches on mode 1

    # Which modes a node may emit has to survive the graph, because downstream a missing transition
    # has two irreconcilable causes -- forbidden by the language, or a successor that escaped the
    # modelled region -- and a universal answer must treat them oppositely.
    @test PCLF.enabled_modes(H) == Dict(1 => Set([1, 2]), 2 => Set([2]))

    no_11 = PCLF.edgeList_to_LabDigraph([(0, 1, 1), (0, 0, 2), (1, 0, 2)])
    @test PCLF.enabled_modes(no_11) == Dict(0 => Set([1, 2]), 1 => Set([2]))

    # A De Bruijn node records history only: every node admits every future word, so a result at one
    # node is a result at all of them. A genuine language restriction is the opposite and the two
    # must be distinguishable without inspecting edges by hand.
    @test !PCLF.restricts_future(de_bruijn(2, 1))
    @test !PCLF.restricts_future(de_bruijn(2, 2))
    @test PCLF.restricts_future(no_11)
    @test PCLF.restricts_future(H)

    for M in (2, 4)
        edges = [(s, mod(s + m - 1, M), m) for s in 0:(M - 1) for m in 1:M]
        R = PCLF.edgeList_to_LabDigraph(edges)
        @test PCLF.is_complete(R, 1:M)
        @test PCLF.is_deterministic(R, 1:M)
    end

    # Co-completeness is the mirror property, and the two separate the De Bruijn graphs exactly:
    # the primal records the last mode played, so every node emits every mode but is entered by only
    # one; the dual commits to the mode to come, so every node emits one mode but is entered by all.
    # Each therefore licenses exactly ONE closed-form common Lyapunov function.
    for k in 1:2
        primal = de_bruijn(2, k)
        dual = de_bruijn(2, k; dual = true)
        @test PCLF.is_complete(primal, 1:2) && !PCLF.is_co_complete(primal, 1:2)
        @test PCLF.is_co_complete(dual, 1:2) && !PCLF.is_complete(dual, 1:2)
    end

    # k = 0 is one node with a self-loop per mode, so it is trivially both.
    single = de_bruijn(3, 0)
    @test PCLF.is_complete(single, 1:3) && PCLF.is_co_complete(single, 1:3)

    # The two properties are independent, not opposites: `H` is entered by both modes at both nodes
    # while node 2 emits only mode 2, so it is co-complete WITHOUT being complete.
    @test PCLF.is_co_complete(H, 1:2) && !PCLF.is_complete(H, 1:2)

    # And a graph can be neither: drop the only mode-1 edge into node 2 and nothing supplies it.
    neither = PCLF.edgeList_to_LabDigraph([(1, 1, 1), (1, 2, 2), (2, 1, 2)])
    @test !PCLF.is_complete(neither, 1:2) && !PCLF.is_co_complete(neither, 1:2)
end

# ------------------------------------------------------------------
# The induced common: which closed form, and no redundant pieces
# ------------------------------------------------------------------
# `build_common_lyapunov` used to always run the observer/subset construction. That is not wrong, but
# on a structured graph it enumerates subsets that cannot contribute: the sublevel set is a UNION over
# observer states of an INTERSECTION over each state's nodes, so any state CONTAINING another is
# redundant. On a complete graph the observer reaches both the full vertex set and every singleton,
# and the full set's intersection then has to be carved out of each singleton -- turning two convex
# polytopes into five disjoint parts, a fragmentation every quotient cell downstream inherits.
@testset "induced common picks the closed form its graph licenses" begin
    mk(G) = PCLF.PCLF(
        G,
        Dict(v => PCLF.PolyhedralPiece([1.0 0.0; 0.0 1.0], [1.0, 1.0]) for v in G.verts),
        0.5,
    )
    states(p) = p.pieces[:clf].observer_states

    primal = de_bruijn(2, 1)
    dual = de_bruijn(2, 1; dual = true)

    # Complete -> V_min: one singleton per node, so the sublevel set is the UNION.
    @test Set(states(PCLF.build_common_lyapunov(mk(primal)))) ==
          Set([Set([(1,)]), Set([(2,)])])
    # Co-complete -> V_max: a single state holding every node, so it is the INTERSECTION.
    @test states(PCLF.build_common_lyapunov(mk(dual))) == [Set([(1,), (2,)])]

    # The observer route must reach the SAME function on the primal graph -- it differs only by
    # enumerating the redundant full set, which the filter then removes.
    @test Set(states(PCLF.build_common_lyapunov(mk(primal); mode = :observer))) ==
          Set([Set([(1,)]), Set([(2,)])])

    # A mode whose structural precondition fails must be refused, not silently returned: min_i V_i is
    # simply not a common Lyapunov function on a graph that is not complete.
    @test_throws ErrorException PCLF.build_common_lyapunov(mk(dual); mode = :min)
    @test_throws ErrorException PCLF.build_common_lyapunov(mk(primal); mode = :max)
    @test_throws ErrorException PCLF.build_common_lyapunov(mk(primal); mode = :nonsense)

    # The filter itself, independent of any graph.
    @test PCLF._drop_redundant_supersets([Set([1, 2]), Set([1]), Set([2])]) ==
          [Set([1]), Set([2])]
    @test PCLF._drop_redundant_supersets([Set([1]), Set([2])]) == [Set([1]), Set([2])]
    @test PCLF._drop_redundant_supersets([Set([1, 2])]) == [Set([1, 2])]
end

@testset "mode_matrices accepts both resetmap shapes" begin
    A1, A2 = [0.5 0.0; 0.0 0.4], [0.3 0.1; 0.0 0.2]
    f = HybridSystems.discreteswitchedsystem([A1, A2])

    # `discreteswitchedsystem` wraps a bare matrix in a `LinearMap`, so this exercises the
    # `:A in fieldnames` branch; the bare-matrix branch is exercised by passing them through.
    @test UT.mode_matrices(f) == [A1, A2]
    @test UT.mode_matrices(f) isa Vector{Matrix{Float64}}
end

end # module TestMain
