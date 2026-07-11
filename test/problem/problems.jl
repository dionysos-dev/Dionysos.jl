module TestMain

using Test
using StaticArrays
using Plots
using Dionysos
import MathematicalSystems

const DI = Dionysos
const UT = DI.Utils
const PR = DI.Problem
const MP = DI.Mapping

sleep(0.1)
println("Started problem tests")

# Minimal dummy system: the problem plot recipes require `MS.stateset(system)`.
struct DummySystem{X}
    X::X
end
MathematicalSystems.stateset(s::DummySystem) = s.X

@testset "Problems" begin
    @testset "Infinity" begin
        @test Base.isfinite(4) == true
        @test Base.isfinite(-2.5) == true
        @test Base.isfinite(Inf) == false

        @test Base.isfinite(PR.Infinity()) == false
    end

    @testset "AlternatingSimulationProblem fields" begin
        X = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))
        sys = DummySystem(X)
        state_set = UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5))

        p = PR.AlternatingSimulationProblem(sys, state_set)

        @test p isa PR.ProblemType
        @test p isa PR.AbstractionProblem
        @test p.system === sys
        @test p.state_set == state_set

        # problems are immutable value objects
        state_set2 = UT.box(SVector(-0.2, -0.2), SVector(0.2, 0.2))
        @test_throws ErrorException p.state_set = state_set2
    end

    @testset "OptimalControlProblem fields" begin
        X = UT.box(SVector(-2.0, -2.0), SVector(2.0, 2.0))
        sys = DummySystem(X)

        XI = UT.box(SVector(-1.0, -1.0), SVector(-0.5, -0.5))
        XT = UT.box(SVector(0.5, 0.5), SVector(1.0, 1.0))

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
        X = UT.box(SVector(-3.0, -3.0), SVector(3.0, 3.0))
        sys = DummySystem(X)

        XI = UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5))
        XS = UT.box(SVector(-2.0, -2.0), SVector(2.0, 2.0))
        T = 7

        p = PR.SafetyProblem(sys, XI, XS, T)

        @test p isa PR.ProblemType
        @test p.system === sys
        @test p.initial_set == XI
        @test p.safe_set == XS
        @test p.time == T
    end

    @testset "CoSafeLTLProblem fields" begin
        X = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))
        sys = DummySystem(X)

        XI = UT.box(SVector(-0.9, -0.9), SVector(-0.8, -0.8))

        # SPEC can be anything; use a Symbol or String as stand-in
        spec = :F_target

        # Label payload can be sets (concrete) or abstract states;
        # here use sets to match your recipe expectation.
        lab = Dict{Symbol, Any}(
            :goal => UT.box(SVector(0.4, 0.4), SVector(0.6, 0.6)),
            :avoid => UT.box(SVector(-0.2, -0.2), SVector(0.2, 0.2)),
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

    @testset "Problem hierarchy and interface" begin
        X = UT.box(SVector(-2.0, -2.0), SVector(2.0, 2.0))
        sys = DummySystem(X)
        XI = UT.box(SVector(-1.0, -1.0), SVector(-0.5, -0.5))
        XT = UT.box(SVector(0.5, 0.5), SVector(1.0, 1.0))

        p = PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 5)
        @test p isa PR.ControlProblem
        @test !(p isa PR.AbstractionProblem)

        # reach-avoid rounds the horizon down; the others round up
        @test PR.horizon_round_up(p) == false
        XS = UT.box(SVector(-1.5, -1.5), SVector(1.5, 1.5))
        @test PR.horizon_round_up(PR.SafetyProblem(sys, XI, XS, 5)) == true

        # remake swaps fields, re-infers type parameters, copies the rest
        XT2 = UT.box(SVector(0.6, 0.6), SVector(1.1, 1.1))
        p2 = PR.remake(p; target_set = XT2, time = 9)
        @test p2 isa PR.OptimalControlProblem
        @test p2.target_set == XT2
        @test p2.time == 9
        @test p2.initial_set === p.initial_set
        @test p2.system === p.system
    end

    @testset "Plots/recipes smoke tests" begin
        # These tests aim for coverage: ensure recipes don't throw.
        # Don't over-assert details to avoid brittleness.

        X = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))
        sys = DummySystem(X)

        region = UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5))
        p_empty = PR.AlternatingSimulationProblem(sys, region)

        XI = UT.box(SVector(-0.9, -0.9), SVector(-0.8, -0.8))
        XT = UT.box(SVector(0.8, 0.8), SVector(0.9, 0.9))
        p_opt = PR.OptimalControlProblem(sys, XI, XT, x -> 0.0, (x, u) -> 1.0, 3)

        XS = UT.box(SVector(-0.7, -0.7), SVector(0.7, 0.7))
        p_safe = PR.SafetyProblem(sys, XI, XS, 4)

        lab = Dict{Symbol, Any}(:ap => UT.box(SVector(0.2, 0.2), SVector(0.3, 0.3)))
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
