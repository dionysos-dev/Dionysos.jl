module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
using IntervalArithmetic
import MathOptInterface as MOI
import MathematicalSystems as MS
import HybridSystems
import LazySets

# Phase 7: where the dynamics come from. They can be written as equations and compiled by a
# backend, or handed over as a Julia function (issue #510).

const WR = Dionysos.Wrapper

function grid_options!(model)
    set_attribute(model, "time_step", 1.0)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    return set_attribute(model, "print_level", 0)
end

@testset "dynamics supplied as a Julia function (#510)" begin
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0, start = -0.75)
    @variable(model, -1.0 <= u <= 1.0)

    # No `∂` anywhere: the states are named instead, the equations are a Julia function, and the
    # time domain is named too — `f(x, u)` alone does not say whether it is a vector field or a
    # one-step map.
    set_role!(x, Dionysos.STATE)
    set_attribute(model, "dynamics", (x, u) -> u)
    set_attribute(model, "time_domain", Dionysos.DISCRETE)
    @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    grid_options!(model)

    optimize!(model)

    @test is_solved_and_feasible(model)
    @test controller_admissible(get_attribute(model, "concrete_controller"), SVector(-0.75))

    # `u` was never declared, so it fell through to being an input; the bounds still come from
    # the model, not from the function.
    problem = get_attribute(model, "concrete_problem")
    @test LazySets.dim(problem.system.X) == 1
    @test LazySets.dim(problem.system.U) == 1
    @test LazySets.low(problem.system.U, 1) == -1.0
end

@testset "set_role! accepts an array and rejects roles with no lowering" begin
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x[1:2] <= 1.0, start = 0.0)
    @variable(model, -1.0 <= u <= 1.0)

    set_role!(x, Dionysos.STATE)        # a whole array at once
    set_attribute(model, "dynamics", (x, u) -> [x[2], u[1]])
    set_attribute(model, "time_domain", Dionysos.DISCRETE)

    opt = backend(model)
    WR.lower(opt)
    @test WR.state_indices(opt.ir) == [1, 2]
    @test WR.input_indices(opt.ir) == [3]

    @test_throws ErrorException set_role!(u, Dionysos.PARAMETER)
end

@testset "a mode can carry its own supplied dynamics" begin
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    @variable(model, -1.0 <= x <= 1.0, start = -0.75)
    @variable(model, -1.0 <= u <= 1.0)
    set_role!(x, Dionysos.STATE)

    @mode(model, slow)
    @mode(model, fast)
    set_attribute(model, "time_domain", Dionysos.DISCRETE)
    set_attribute(slow, "dynamics", (x, u) -> 0.5 .* u)
    set_attribute(fast, "dynamics", (x, u) -> 2.0 .* u)
    @constraint(fast, [x] in Final(UT.box(SVector(-0.5), SVector(0.5))))

    add_transition!(model, slow => fast) do t
        return @constraint(t, x >= 0.0)
    end

    opt = backend(model)
    problem = WR.lower(opt)

    hs = problem.system
    @test MS.mapping(HybridSystems.mode(hs, 1))([0.0], [1.0])[1] ≈ 0.5
    @test MS.mapping(HybridSystems.mode(hs, 2))([0.0], [1.0])[1] ≈ 2.0
end

@testset "an explicitly declared clock wins over the inference" begin
    # With supplied dynamics there are no equations to recognise a clock from, so it is named.
    # The inference must then stand down rather than labelling a second variable as the clock.
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    @variable(model, -1.0 <= x <= 1.0, start = 0.0)
    @variable(model, -1.0 <= u <= 1.0)
    @variable(model, 0.0 <= t <= 5.0)

    set_role!(x, Dionysos.STATE)
    set_role!(t, Dionysos.CLOCK)

    @mode(model, a)
    @mode(model, b)
    set_attribute(a, "dynamics", (x, u) -> u)
    set_attribute(b, "dynamics", (x, u) -> u)
    @constraint(a, ∂(t) == 1)
    @constraint(b, ∂(t) == 1)
    @constraint(b, [x] in Final(UT.box(SVector(-0.5), SVector(0.5))))
    add_transition!(model, a => b) do tr
        return @constraint(tr, x >= 0.0)
    end

    opt = backend(model)
    WR.lower(opt)

    @test WR.clock_index(opt.ir) == 3
    @test count(v -> v.role === Dionysos.CLOCK, opt.ir.variables) == 1
    @test WR.state_indices(opt.ir) == [1]
end

@testset "the evaluator backend needs no optional dependency" begin
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0, start = -0.75)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == u)
    @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    grid_options!(model)
    set_attribute(model, "dynamics_backend", WR.NonlinearEvaluatorBackend())

    optimize!(model)

    @test is_solved_and_feasible(model)
    @test controller_admissible(get_attribute(model, "concrete_controller"), SVector(-0.75))
end

@testset "the evaluator and symbolic backends agree" begin
    function compiled(backend)
        model = direct_model(Dionysos.Optimizer())
        @variable(model, -2.0 <= x[1:2] <= 2.0)
        @variable(model, -1.0 <= u <= 1.0)
        @constraint(model, ∂(x[1]) == x[2] * cos(x[1]))
        @constraint(model, ∂(x[2]) == u - 0.5 * x[2])
        opt = backend_of(model)
        WR._apply_options!(opt)
        WR.infer_roles!(opt.ir)
        return WR.compile_dynamics(backend, opt.ir, opt.ir.dynamics)
    end
    backend_of(m) = JuMP.backend(m)

    symbolic = compiled(WR.SymbolicADBackend())
    evaluated = compiled(WR.NonlinearEvaluatorBackend())

    for (x, u) in (
        (SVector(0.3, -0.7), SVector(0.2)),
        (SVector(-1.1, 1.4), SVector(-0.9)),
        (SVector(0.0, 0.0), SVector(0.0)),
    )
        @test collect(symbolic(x, u)) ≈ collect(evaluated(x, u))
    end
end

@testset "the growth bound is derived when none is given" begin
    # `System.compute_jacobian_bound` traces the dynamics and bounds the Jacobian over X with
    # interval arithmetic, so `GROWTH` does not need a hand-written `jacobian_bound`.
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0, start = -0.75)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == -0.5 * x + u)
    @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))

    set_attribute(model, "time_step", 0.5)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "print_level", 0)
    # Deliberately no `jacobian_bound`.

    optimize!(model)
    @test is_solved_and_feasible(model)
end

# `f(x, u)` is equally readable as a vector field or as a one-step map, and the two describe
# different plants. Before this was checked, a supplied-dynamics model lowered silently as
# discrete-time — it built and simulated, just not the system the user wrote.
@testset "supplied dynamics must declare their time domain" begin
    function supplied_model(; domain = nothing)
        model = direct_model(Dionysos.Optimizer())
        @variable(model, -1.0 <= x <= 1.0, start = -0.75)
        @variable(model, -1.0 <= u <= 1.0)
        set_role!(x, Dionysos.STATE)
        set_attribute(model, "dynamics", (x, u) -> u)
        domain === nothing || set_attribute(model, "time_domain", domain)
        @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
        return model
    end

    # Unknown: refused, and the message says what to write.
    err = try
        WR.lower(backend(supplied_model()))
        nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing
    @test occursin("time domain is unknown", err)
    @test occursin("time_domain", err)

    @test WR.lower(backend(supplied_model(; domain = Dionysos.CONTINUOUS))).system isa
          MS.ConstrainedBlackBoxControlContinuousSystem
    @test WR.lower(backend(supplied_model(; domain = Dionysos.DISCRETE))).system isa
          MS.ConstrainedBlackBoxControlDiscreteSystem

    # Written dynamics already fix the domain; an attribute contradicting them is refused
    # rather than silently overriding what the model says.
    written = direct_model(Dionysos.Optimizer())
    @variable(written, -1.0 <= y <= 1.0, start = 0.0)
    @variable(written, -1.0 <= v <= 1.0)
    @constraint(written, ∂(y) == v)
    set_attribute(written, "time_domain", Dionysos.DISCRETE)
    @constraint(written, final(y) in MOI.Interval(-0.5, 0.5))
    @test_throws ErrorException WR.lower(backend(written))
end

end # module TestMain
