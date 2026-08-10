module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
import MathOptInterface as MOI
import MathematicalSystems as MS
import HybridSystems
import LazySets

# Obstacles and guards beyond boxes: `∉` with any bounded `LazySet`, guards written as a
# `Guard(S)` set or as a multi-variable affine constraint, and the cases that stay rejected.

const WR = Dionysos.Wrapper

function integrator_2d!(model)
    @variable(model, -2.0 <= x <= 2.0, start = -1.5)
    @variable(model, -2.0 <= y <= 2.0, start = -1.5)
    @variable(model, -1.0 <= u[1:2] <= 1.0)
    @constraint(model, ∂(x) == u[1])
    @constraint(model, ∂(y) == u[2])
    return model
end

function lowered_problem(build!)
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    build!(model)
    return WR.lower(backend(model))
end

# ----------------------------------------------------------------------------------------
# Obstacles
# ----------------------------------------------------------------------------------------

@testset "an obstacle may be any bounded LazySet" begin
    ball = LazySets.Ball2([0.0, 0.0], 0.5)
    problem = lowered_problem() do model
        integrator_2d!(model)
        @constraint(model, [model[:x], model[:y]] ∉ ball)
        return @constraint(
            model,
            [model[:x], model[:y]] in Final(
                LazySets.Hyperrectangle(;
                    low = SVector(1.0, 1.0),
                    high = SVector(1.5, 1.5),
                ),
            )
        )
    end

    X = problem.system.X
    @test SVector(0.0, 0.0) ∉ X          # carved out by the ball
    @test SVector(1.2, 1.2) ∈ X
    @test SVector(1.9, 1.9) ∈ X
end

@testset "a box obstacle still spans the coordinates it omits" begin
    # `MOI.HyperRectangle` keeps the old behaviour: written over a subset of the coordinates,
    # extruded across the variable bounds on the rest.
    problem = lowered_problem() do model
        integrator_2d!(model)
        @constraint(model, [model[:x]] ∉ MOI.HyperRectangle([-0.25], [0.25]))
        return @constraint(
            model,
            [model[:x], model[:y]] in Final(
                LazySets.Hyperrectangle(;
                    low = SVector(1.0, 1.0),
                    high = SVector(1.5, 1.5),
                ),
            )
        )
    end

    X = problem.system.X
    @test SVector(0.0, -1.9) ∉ X         # the strip spans all of `y`
    @test SVector(0.0, 1.9) ∉ X
    @test SVector(1.0, 0.0) ∈ X
end

@testset "a non-box obstacle must span the whole state vector" begin
    ball = LazySets.Ball2([0.0], 0.5)
    err = try
        lowered_problem() do model
            integrator_2d!(model)
            @constraint(model, [model[:x]] ∉ ball)
            return @constraint(
                model,
                [model[:x], model[:y]] in Final(
                    LazySets.Hyperrectangle(;
                        low = SVector(1.0, 1.0),
                        high = SVector(1.5, 1.5),
                    ),
                )
            )
        end
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("whole state vector", sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------
# Guards
# ----------------------------------------------------------------------------------------

# Two modes with the same dynamics, so only the guard distinguishes the models below.
function two_mode_model(guard!)
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    integrator_2d!(model)

    @mode(model, a)
    @mode(model, b)
    for m in (a, b)
        @constraint(m, ∂(model[:x]) == model[:u][1])
        @constraint(m, ∂(model[:y]) == model[:u][2])
    end
    @constraint(
        b,
        [model[:x], model[:y]] in
        Final(LazySets.Hyperrectangle(; low = SVector(1.0, 1.0), high = SVector(1.5, 1.5)))
    )

    add_transition!(model, a => b) do t
        return guard!(t, model)
    end
    return model
end

function _guard_of(problem)
    hs = problem.system
    transition = first(HybridSystems.transitions(hs.automaton))
    return MS.stateset(HybridSystems.resetmap(hs, transition))
end

@testset "a guard may be a LazySet" begin
    ball = LazySets.Ball2([0.0, 0.0], 1.0)
    model = two_mode_model() do t, m
        return @constraint(t, [m[:x], m[:y]] in Guard(ball))
    end
    guard = _guard_of(WR.lower(backend(model)))

    @test SVector(0.0, 0.0) ∈ guard
    @test SVector(0.5, 0.5) ∈ guard
    @test SVector(1.9, 1.9) ∉ guard      # outside the ball, inside the mode's box
end

@testset "a guard may be a multi-variable affine constraint" begin
    model = two_mode_model() do t, m
        return @constraint(t, m[:x] + m[:y] <= 1.0)
    end
    guard = _guard_of(WR.lower(backend(model)))

    @test SVector(0.0, 0.0) ∈ guard
    @test SVector(-1.0, 1.0) ∈ guard     # on the boundary of the half-space
    @test SVector(1.5, 1.5) ∉ guard
end

@testset "guards intersect: bounds, half-spaces and sets together" begin
    ball = LazySets.Ball2([0.0, 0.0], 1.5)
    model = two_mode_model() do t, m
        @constraint(t, m[:x] >= 0.0)                       # per-coordinate bound
        @constraint(t, m[:x] + m[:y] <= 1.0)               # half-space
        return @constraint(t, [m[:x], m[:y]] in Guard(ball))  # set
    end
    guard = _guard_of(WR.lower(backend(model)))

    @test SVector(0.25, 0.25) ∈ guard    # satisfies all three
    @test SVector(-0.5, 0.0) ∉ guard     # fails the bound
    @test SVector(1.0, 1.0) ∉ guard      # fails the half-space
    @test SVector(0.1, -1.9) ∉ guard     # fails the ball
end

@testset "a plain box guard is still a plain box" begin
    # The fast path must not change: with only per-coordinate bounds the guard stays a
    # `Hyperrectangle`, not a lazy intersection.
    model = two_mode_model() do t, m
        return @constraint(t, m[:x] >= 1.0)
    end
    guard = _guard_of(WR.lower(backend(model)))

    @test guard isa LazySets.AbstractHyperrectangle
    @test LazySets.low(guard, 1) == 1.0
end

@testset "a Guard set must span the state vector, and belongs on a transition" begin
    small = LazySets.Ball2([0.0], 1.0)
    model = two_mode_model() do t, m
        return @constraint(t, [m[:x]] in Guard(small))
    end
    err = try
        WR.lower(backend(model))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("whole state vector", sprint(showerror, err))

    # On a mode, `Guard` is the wrong marker.
    model2 = direct_model(Dionysos.Optimizer())
    integrator_2d!(model2)
    @mode(model2, only_mode)
    @test_throws ErrorException @constraint(
        only_mode,
        [model2[:x], model2[:y]] in Guard(LazySets.Ball2([0.0, 0.0], 1.0))
    )
end

@testset "a multi-variable constraint on a mode is still rejected" begin
    model = direct_model(Dionysos.Optimizer())
    integrator_2d!(model)
    @mode(model, m1)
    @test_throws ErrorException @constraint(m1, model[:x] + model[:y] <= 1.0)
end

@testset "a non-box guard is rejected on a clocked model" begin
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    integrator_2d!(model)
    @variable(model, 0.0 <= clk <= 5.0, start = 0.0)
    Dionysos.set_role!(clk, Dionysos.CLOCK)

    @mode(model, a)
    @mode(model, b)
    for m in (a, b)
        @constraint(m, ∂(model[:x]) == model[:u][1])
        @constraint(m, ∂(model[:y]) == model[:u][2])
        @constraint(m, ∂(clk) == 1)
    end
    @constraint(
        b,
        [model[:x], model[:y]] in
        Final(LazySets.Hyperrectangle(; low = SVector(1.0, 1.0), high = SVector(1.5, 1.5)))
    )
    add_transition!(model, a => b) do t
        return @constraint(t, model[:x] + model[:y] <= 1.0)
    end

    err = try
        WR.lower(backend(model))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("clock", sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------
# Switching and objective
# ----------------------------------------------------------------------------------------

@testset "a degenerate guard is diagnosed, not asserted from inside LazySets" begin
    # Coefficients that cancel leave no state in the constraint, so there is no half-space to
    # build. `x - x <= 1` is vacuous; `x - x <= -1` can never hold.
    vacuous = two_mode_model() do t, m
        @constraint(t, m[:x] - m[:x] <= 1.0)
        return @constraint(t, m[:x] >= 0.0)        # a real guard alongside it
    end
    guard = _guard_of(WR.lower(backend(vacuous)))
    @test SVector(0.5, 0.0) ∈ guard                # the vacuous constraint changed nothing
    @test SVector(-0.5, 0.0) ∉ guard

    unsatisfiable = two_mode_model() do t, m
        return @constraint(t, m[:x] - m[:x] <= -1.0)
    end
    err = try
        WR.lower(backend(unsatisfiable))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    message = sprint(showerror, err)
    @test occursin("cancel out", message)
    @test occursin("Transition", message)          # names the transition, unlike the assertion
end

@testset "contradictory guard bounds are diagnosed" begin
    # `Hyperrectangle` would build a negative radius and assert from inside LazySets, naming neither
    # the transition nor the variable.
    model = two_mode_model() do t, m
        @constraint(t, m[:x] >= 1.0)
        return @constraint(t, m[:x] <= 0.0)
    end
    err = try
        WR.lower(backend(model))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    message = sprint(showerror, err)
    @test occursin("empty guard", message)
    @test occursin("x", message)                   # names the offending variable

    # A guard entirely outside the source mode's own bounds is empty for the same reason, and
    # that intersection is the non-obvious part, so the message must mention it.
    outside = two_mode_model() do t, m
        return @constraint(t, m[:x] >= 99.0)
    end
    outside_err = try
        WR.lower(backend(outside))
        nothing
    catch e
        e
    end
    @test outside_err isa ErrorException
    @test occursin("source mode", sprint(showerror, outside_err))
end

@testset "an unsupported obstacle set is refused clearly" begin
    # `∉ MOI.Nonnegatives(2)` used to reach the lowering and raise a bare `MethodError` from
    # `LazySets.dim`. JuMP now turns it away as an unsupported constraint.
    model = direct_model(Dionysos.Optimizer())
    integrator_2d!(model)
    @test_throws Exception @constraint(model, [model[:x], model[:y]] ∉ MOI.Nonnegatives(2))
end

# The end-to-end counterpart — solving a model whose guard is a half-space — lives in
# `hybrid.jl`, which is tagged `:slow`. Everything above inspects the lowered problem instead,
# so this file stays cheap enough for the fast suite.

@testset "controlled switching is refused, not silently ignored" begin
    model = direct_model(Dionysos.Optimizer())
    integrator_2d!(model)
    @mode(model, a)
    @mode(model, b)

    err = try
        add_transition!(model, a => b; switching = HybridSystems.ControlledSwitching()) do t
            return @constraint(t, model[:x] >= 1.0)
        end
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("AutonomousSwitching", sprint(showerror, err))

    # The default is autonomous and still works.
    @test add_transition!(model, a => b) do t
        return @constraint(t, model[:x] >= 1.0)
    end isa WR.Transition
end

@testset "the @objective error points somewhere that exists" begin
    err = try
        lowered_problem() do model
            integrator_2d!(model)
            @constraint(
                model,
                [model[:x], model[:y]] in Final(
                    LazySets.Hyperrectangle(;
                        low = SVector(1.0, 1.0),
                        high = SVector(1.5, 1.5),
                    ),
                )
            )
            return @objective(model, Min, model[:u][1])
        end
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    message = sprint(showerror, err)
    # It used to name `set_attribute(model, "transition_cost", …)`, which throws
    # `MOI.UnsupportedAttribute`: the optimizer has no such field.
    @test !occursin("set_attribute", message)
    @test occursin("concrete_problem", message)
end

end # module TestMain
