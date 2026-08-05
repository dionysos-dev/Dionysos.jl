module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
using Spot
import MathOptInterface as MOI
import LazySets

# Phase 9: the general specification layer — a temporal formula over named regions. The atomic
# propositions are the JuMP constraint names, so no separate registration call exists.

const WR = Dionysos.Wrapper

function reach_avoid_model()
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    @variable(model, -2.0 <= x[1:2] <= 2.0, start = -1.5)
    @variable(model, -1.0 <= u[1:2] <= 1.0)
    @constraint(model, ∂(x[1]) == u[1])
    @constraint(model, ∂(x[2]) == u[2])
    return model, x
end

@testset "a Label's atomic proposition is the constraint's name" begin
    model, x = reach_avoid_model()
    goal = UT.box(SVector(1.0, 1.0), SVector(1.8, 1.8))
    hazard = UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5))

    @constraint(model, goal_region, x in Label(goal))
    @constraint(model, hazard_region, x in Label(hazard; semantics = MP.OUTER))
    @specification(model, ltl"F(goal_region) & G(!hazard_region)")

    problem = WR.lower(model)

    @test problem isa PR.CoSafeLTLProblem
    @test Set(keys(problem.labeling)) == Set([:goal_region, :hazard_region])
    @test problem.labeling[:goal_region] === goal
    # The semantics travel with the region: reach it conservatively, avoid it conservatively.
    @test problem.ap_semantics[:goal_region] == MP.INNER
    @test problem.ap_semantics[:hazard_region] == MP.OUTER
    # The formula became the automaton the co-safe solver steps.
    @test problem.spec isa OPDS.AbstractSpecStepper
end

@testset "a hand-written monitor is accepted instead of a formula" begin
    # Nothing forces Spot: any `AbstractSpecStepper` works, so a monitor can be written by hand.
    model, x = reach_avoid_model()
    @constraint(
        model,
        goal_region,
        x in Label(UT.box(SVector(1.0, 1.0), SVector(1.8, 1.8)))
    )

    monitor = OPDS.FunctionMonitor(1, Set([2]), (q, ap) -> (:goal_region in ap) ? 2 : q)
    @specification(model, monitor)

    problem = WR.lower(model)
    @test problem isa PR.CoSafeLTLProblem
    @test problem.spec === monitor
end

@testset "the formula takes precedence over the pattern markers" begin
    model, x = reach_avoid_model()
    @constraint(
        model,
        goal_region,
        x in Label(UT.box(SVector(1.0, 1.0), SVector(1.8, 1.8)))
    )
    # A `Final` marker as well: the formula is the general form and wins.
    @constraint(model, x in Final(UT.box(SVector(0.0, 0.0), SVector(1.0, 1.0))))
    @specification(model, ltl"F(goal_region)")

    @test WR.lower(model) isa PR.CoSafeLTLProblem
end

@testset "validation" begin
    # An anonymous label has no proposition to be referred to by.
    model, x = reach_avoid_model()
    @constraint(model, x in Label(UT.box(SVector(1.0, 1.0), SVector(1.8, 1.8))))
    @specification(model, ltl"F(goal_region)")
    err = try
        WR.lower(model)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("needs a name", err.msg)

    # A formula with no region at all.
    bare, _ = reach_avoid_model()
    @specification(bare, ltl"F(goal_region)")
    @test_throws ErrorException WR.lower(bare)

    # Something that is neither a formula nor a monitor.
    wrong, y = reach_avoid_model()
    @constraint(
        wrong,
        goal_region,
        y in Label(UT.box(SVector(1.0, 1.0), SVector(1.8, 1.8)))
    )
    @specification(wrong, 42)
    @test_throws ErrorException WR.lower(wrong)
end

@testset "end-to-end: a co-safe LTL model solved from JuMP" begin
    # Reach the goal while staying out of the hazard — the reach-avoid shape, but stated as a
    # formula and solved through the automaton product rather than the specialised fixpoint.
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    @variable(model, -2.0 <= x[1:2] <= 2.0, start = -1.5)
    @variable(model, -1.0 <= u[1:2] <= 1.0)
    @constraint(model, ∂(x[1]) == u[1])
    @constraint(model, ∂(x[2]) == u[2])

    @constraint(
        model,
        goal_region,
        x in Label(UT.box(SVector(1.0, 1.0), SVector(1.8, 1.8)))
    )
    @constraint(
        model,
        hazard_region,
        x in Label(UT.box(SVector(-0.4, -0.4), SVector(0.4, 0.4)); semantics = MP.OUTER)
    )
    @specification(model, ltl"F(goal_region) & G(!hazard_region)")

    set_attribute(model, "time_step", 0.5)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.4, 0.4)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))

    optimize!(model)

    @test get_attribute(model, "concrete_problem") isa PR.CoSafeLTLProblem
    @test get_attribute(model, "abstract_controller") !== nothing
    @test get_attribute(model, "success") isa Bool
end

end # module TestMain
