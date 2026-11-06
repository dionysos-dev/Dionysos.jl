module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import PathCompleteCertificates as PCC
import Clarabel
import HiGHS
import JuMP
import LazySets
import LinearAlgebra as LA
import Logging
import MathOptInterface as MOI
import Random

const PCLF = UT.PathCompleteFramework
const PCQ = AB.PCLFBisimulationQuotient

# A stable pair, certified by every template below.
const A = [[0.70 0.10; 0.00 0.65], [0.60 -0.15; 0.10 0.55]]
const PROBLEM = PCC.StabilityProblem(PCC.switched_system(A))
const LP = JuMP.optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
const SDP = JuMP.optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

const SAMPLES = let rng = Random.MersenneTwister(1)
    [randn(rng, 2) for _ in 1:200]
end

# The converted pieces must be the certificate's node functions, not approximations of them.
# `symmetrized` compares with max(V(x), V(-x)), what a symmetric Dionysos piece can hold.
function same_functions(pclf, certificate; symmetrized = false)
    value(V, x) = PCC.node_value(PROBLEM, V, x)
    expected(V, x) = symmetrized ? max(value(V, x), value(V, -x)) : value(V, x)
    return all(
        isapprox(PCLF.piece_value(pclf.pieces[k], x), expected(V, x); rtol = 1e-8) for
        (k, V) in enumerate(PCC.functions(certificate)) for x in SAMPLES
    )
end

# Every edge inequality V_d(A_m x) ≤ γ^δ V_s(x) of the converted certificate, at its own rate.
function is_certificate(pclf, δ)
    γ = pclf.JSRapprox
    return all(
        PCLF.piece_value(pclf.pieces[d], A[m] * x) <=
        γ^δ * PCLF.piece_value(pclf.pieces[s], x) * (1 + 1e-6) for
        (s, d, m) in pclf.graph.edges for x in SAMPLES
    )
end

@testset "dual piecewise linear: faces normal to the rows of G" begin
    θ = π / 6
    Rot = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    certificate = PCC.certify(
        PCC.DualPiecewiseLinearTemplate([Rot, Rot]),
        PCC.de_bruijn(1, 2),
        PROBLEM;
        optimizer = LP,
    )
    pclf = PCLF.PCLF(certificate)

    @test pclf.JSRapprox == certificate.rate
    @test pclf.graph.verts == Set([1, 2])
    @test Set(pclf.graph.edges) == Set([(1, 1, 1), (1, 2, 2), (2, 1, 1), (2, 2, 2)])
    @test same_functions(pclf, certificate)
    @test is_certificate(pclf, 1)
    for piece in values(pclf.pieces)
        @test piece isa PCLF.PolyhedralPiece
        @test size(piece.G) == (2, 2)
        # Each row is a row of the rotation, up to sign: the faces keep their directions.
        @test all(
            any(abs(LA.dot(g, r)) ≈ 1 for r in eachrow(Rot)) for g in eachrow(piece.G)
        )
        @test all(>(0), piece.w)
        # Two pairs of parallel faces, so exactly four halfspaces.
        @test length(LazySets.constraints_list(PCLF.get_sublevel_set(piece, 1.0))) == 4
    end
end

@testset "piecewise linear: symmetrized, still a certificate" begin
    partition = PCC.planar_conic_partition(2)   # 16 cones
    certificate = PCC.certify(
        PCC.PiecewiseLinearTemplate([partition, partition]),
        PCC.dual(PCC.de_bruijn(1, 2)),
        PROBLEM;
        optimizer = LP,
    )
    pclf = PCLF.PCLF(certificate)

    # The template's functions need not be symmetric, and with HiGHS here they are not.
    @test same_functions(pclf, certificate; symmetrized = true)
    @test is_certificate(pclf, 1)
    for piece in values(pclf.pieces)
        # One row per linear piece, up to sign and repetition.
        @test 2 <= size(piece.G, 1) <= 16
        @test piece.w == ones(size(piece.G, 1))
        @test allunique(eachrow(piece.G))
    end
    # The dual De Bruijn graph is co-complete, so the induced common is the maximum of the pieces.
    common = PCLF.build_common_lyapunov(pclf)
    @test common.pieces[:clf].observer_states == [Set([1, 2])]
end

@testset "quadratic: the matrix of the quadratic form" begin
    certificate =
        PCC.certify(PCC.QuadraticTemplate(), PCC.common_graph(2), PROBLEM; optimizer = SDP)
    pclf = PCLF.PCLF(certificate)

    piece = pclf.pieces[1]
    @test piece isa PCLF.EllipsoidalPiece
    @test piece.P ≈ only(PCC.functions(certificate)).Q
    @test same_functions(pclf, certificate)
    @test is_certificate(pclf, 2)
end

@testset "labels, and the certificate as an optimizer attribute" begin
    certificate = PCC.certify(
        PCC.DualPiecewiseLinearTemplate(2, 2),
        PCC.de_bruijn(1, 2),
        PROBLEM;
        optimizer = LP,
    )
    labelled = PCLF.PCLF(certificate; labels = [(1,), (2,)])
    @test labelled.graph.verts == Set([(1,), (2,)])
    @test Set(keys(labelled.pieces)) == Set([(1,), (2,)])
    @test ((2,), (1,), 1) in labelled.graph.edges

    optimizer = MOI.instantiate(PCQ.OptimizerBisimulationQuotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), certificate)
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("pclf")) isa PCLF.PCLF
end

@testset "what has no geometry is refused" begin
    # A rate below the joint spectral radius has no certificate.
    infeasible = PCC.certify(
        PCC.DualPiecewiseLinearTemplate(1, 2),
        PCC.common_graph(2),
        PROBLEM;
        optimizer = LP,
        rate = 0.1,
    )
    @test !PCC.is_feasible(infeasible)
    @test_throws ArgumentError PCLF.PCLF(infeasible)

    # A copositive node function is defined on the orthant only: no polytope to read. Its program
    # makes Polyhedra warn, many times, that it detects linearity without a solver, which is
    # irrelevant to the conversion.
    positive =
        PCC.StabilityProblem(PCC.switched_system([[0.5 0.1; 0.2 0.4], [0.3 0.2; 0.1 0.5]]))
    copositive = Logging.with_logger(Logging.NullLogger()) do
        return PCC.certify(
            PCC.LinearCopositiveTemplate(),
            PCC.common_graph(2),
            positive;
            optimizer = LP,
        )
    end
    @test_throws ArgumentError PCLF.PCLF(copositive)
end

end # module TestMain
