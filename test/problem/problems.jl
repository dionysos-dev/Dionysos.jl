module TestMain

using Test
using StaticArrays
using Plots
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const PR = DI.Problem
const MP = DI.Mapping

sleep(0.1)
println("Started problem tests")

# Minimal dummy system type: only needs field `X` for recipes
struct DummySystem{X}
    X::X
end

@testset "Problems" begin
    @testset "Infinity" begin
        @test Base.isfinite(4) == true
        @test Base.isfinite(-2.5) == true
        @test Base.isfinite(Inf) == false

        @test Base.isfinite(PR.Infinity()) == false
        # Optional sanity: should be a Real subtype
        @test PR.Infinity() isa Real
    end

    @testset "AlternatingSimulationProblem fields" begin
        X = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))
        sys = DummySystem(X)
        region = UT.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))

        p = PR.AlternatingSimulationProblem(sys, region)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.region == region

        # mutation (it's mutable)
        region2 = UT.HyperRectangle(SVector(-0.2, -0.2), SVector(0.2, 0.2))
        p.region = region2
        @test p.region == region2
    end

    @testset "OptimalControlProblem fields" begin
        X = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
        sys = DummySystem(X)

        XI = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(-0.5, -0.5))
        XT = UT.HyperRectangle(SVector(0.5, 0.5), SVector(1.0, 1.0))

        state_cost = x -> 0.0
        transition_cost = (x, u) -> 1.0
        T = 5

        p = PR.OptimalControlProblem(sys, XI, XT, state_cost, transition_cost, T)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.initial_set == XI
        @test p.target_set == XT
        @test p.state_cost === state_cost
        @test p.transition_cost === transition_cost
        @test p.time == T
    end

    @testset "SafetyProblem fields" begin
        X = UT.HyperRectangle(SVector(-3.0, -3.0), SVector(3.0, 3.0))
        sys = DummySystem(X)

        XI = UT.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))
        XS = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
        T = 7

        p = PR.SafetyProblem(sys, XI, XS, T)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.initial_set == XI
        @test p.safe_set == XS
        @test p.time == T
    end

    @testset "CoSafeLTLProblem fields" begin
        X = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))
        sys = DummySystem(X)

        XI = UT.HyperRectangle(SVector(-0.9, -0.9), SVector(-0.8, -0.8))

        # SPEC can be anything; use a Symbol or String as stand-in
        spec = :F_target

        # Label payload can be sets (concrete) or abstract states;
        # here use sets to match your recipe expectation.
        lab = Dict{Symbol, Any}(
            :goal => UT.HyperRectangle(SVector(0.4, 0.4), SVector(0.6, 0.6)),
            :avoid => UT.HyperRectangle(SVector(-0.2, -0.2), SVector(0.2, 0.2)),
        )

        ap_sem = Dict{Symbol, Any}(:goal => MP.INNER, :avoid => MP.OUTER)

        p = PR.CoSafeLTLProblem(sys, XI, spec, lab, ap_sem)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.initial_set == XI
        @test p.spec == spec
        @test p.labeling == lab
        @test p.ap_semantics == ap_sem
    end

    @testset "Plots/recipes smoke tests" begin
        # These tests aim for coverage: ensure recipes don't throw.
        # Don't over-assert details to avoid brittleness.

        X = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))
        sys = DummySystem(X)

        region = UT.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))
        p_empty = PR.AlternatingSimulationProblem(sys, region)

        XI = UT.HyperRectangle(SVector(-0.9, -0.9), SVector(-0.8, -0.8))
        XT = UT.HyperRectangle(SVector(0.8, 0.8), SVector(0.9, 0.9))
        p_opt = PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 3)

        XS = UT.HyperRectangle(SVector(-0.7, -0.7), SVector(0.7, 0.7))
        p_safe = PR.SafetyProblem(sys, XI, XS, 4)

        lab = Dict{Symbol, Any}(
            :ap => UT.HyperRectangle(SVector(0.2, 0.2), SVector(0.3, 0.3)),
        )
        ap_sem = Dict{Symbol, Any}(:ap => MP.INNER)
        p_ltl = PR.CoSafeLTLProblem(sys, XI, :spec, lab, ap_sem)

        @test_nowarn plot(p_empty)
        @test_nowarn plot(p_opt)
        @test_nowarn plot(p_safe)
        @test_nowarn plot(p_ltl)
    end
end

println("End problem tests")

end # module TestMain
