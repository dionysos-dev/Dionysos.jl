module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
import MathOptInterface as MOI
import LazySets

# Phase 3: the specification markers, the problem type inferred from them, the horizon, and
# the two validation errors (unbounded variables, `@objective`).

const WR = Dionysos.Wrapper

function lowered_problem(build!; time_step = nothing)
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    build!(model)
    opt = backend(model)
    WR._apply_options!(opt)
    WR.infer_roles!(opt.ir)
    f = WR.compile_dynamics(opt.dynamics_backend, opt.ir)
    return WR.build_problem(opt.ir, f; time_step = time_step)
end

# A 1-D single integrator; every test below differs only in its specification. The variables
# are reached through the model's object dictionary (`model[:x]`), since `@variable` binds the
# Julia name inside this function only.
function integrator!(model)
    @variable(model, -2.0 <= x <= 2.0, start = -1.5)
    @variable(model, -1.0 <= u <= 1.0)
    return @constraint(model, ∂(x) == u)
end

state(model) = [model[:x]]

@testset "problem type is inferred from the specification markers" begin
    target = UT.box(SVector(-0.5), SVector(0.5))
    safe = UT.box(SVector(-1.8), SVector(1.8))

    reach = lowered_problem() do model
        integrator!(model)
        return @constraint(model, state(model) in Final(target))
    end
    @test reach isa PR.OptimalControlProblem
    @test reach.target_set === target

    safety = lowered_problem() do model
        integrator!(model)
        return @constraint(model, state(model) in Always(safe))
    end
    @test safety isa PR.SafetyProblem
    @test safety.safe_set === safe

    stay = lowered_problem() do model
        integrator!(model)
        return @constraint(model, state(model) in EventuallyAlways(target))
    end
    @test stay isa PR.ReachAndStayProblem
    @test stay.target_set === target
    # With no `Always`, the safe set defaults to the whole state box.
    @test LazySets.low(stay.safe_set, 1) == -2.0

    stay_safe = lowered_problem() do model
        integrator!(model)
        @constraint(model, state(model) in EventuallyAlways(target))
        return @constraint(model, state(model) in Always(safe))
    end
    @test stay_safe isa PR.ReachAndStayProblem
    @test stay_safe.safe_set === safe

    # No specification at all: the task is the abstraction itself. This was impossible to
    # express from JuMP before.
    abstraction = lowered_problem(integrator!)
    @test abstraction isa PR.AlternatingSimulationProblem
end

@testset "reach-avoid carries the Always set as the problem's safe set" begin
    # The `Always` set travels as `safe_set` rather than being folded into `X`. That keeps the
    # unsafe region representable, so the synthesis can reason about it instead of it simply
    # not existing — the distinction between `Always` and `∉`.
    safe = UT.box(SVector(-1.0), SVector(1.0))
    problem = lowered_problem() do model
        integrator!(model)
        @constraint(model, state(model) in Final(UT.box(SVector(-0.5), SVector(0.5))))
        return @constraint(model, state(model) in Always(safe))
    end

    @test problem isa PR.OptimalControlProblem
    @test problem.safe_set == safe
    # The state space is the declared variable bounds, untouched by the `Always` set.
    @test SVector(0.0) ∈ problem.system.X
    @test SVector(1.5) ∈ problem.system.X
end

@testset "a reach-avoid safe set may be any bounded LazySet" begin
    # Folding into `X` demanded a box. Now that the set has a field of its own, the only
    # requirement is that the discretisation can handle it — which it always could.
    ball = LazySets.Ball2(zeros(1), 1.0)
    problem = lowered_problem() do model
        integrator!(model)
        @constraint(model, state(model) in Final(UT.box(SVector(-0.5), SVector(0.5))))
        return @constraint(model, state(model) in Always(ball))
    end

    @test problem isa PR.OptimalControlProblem
    @test problem.safe_set === ball
end

@testset "a reach model without Always has no safe set" begin
    problem = lowered_problem() do model
        integrator!(model)
        return @constraint(
            model,
            state(model) in Final(UT.box(SVector(-0.5), SVector(0.5)))
        )
    end
    @test problem.safe_set === nothing
end

@testset "the initial set comes from Start, start(...) or the start keyword" begin
    # The `start = v` keyword gives a singleton.
    singleton = lowered_problem() do model
        integrator!(model)
        return @constraint(
            model,
            state(model) in Final(UT.box(SVector(-0.5), SVector(0.5)))
        )
    end
    @test LazySets.low(singleton.initial_set, 1) == -1.5
    @test LazySets.high(singleton.initial_set, 1) == -1.5

    # A `Start` marker gives a region and takes precedence.
    region = UT.box(SVector(-1.9), SVector(-1.1))
    problem = lowered_problem() do model
        integrator!(model)
        @constraint(model, state(model) in Final(UT.box(SVector(-0.5), SVector(0.5))))
        return @constraint(model, state(model) in Start(region))
    end
    @test problem.initial_set === region
end

@testset "specification sets accept any bounded LazySet" begin
    # A ball, not a box — the discretization layer has always accepted these; only the
    # front-end restricted the user to boxes.
    ball = LazySets.Ball2([0.0, 0.0], 0.5)
    problem = lowered_problem() do model
        @variable(model, -2.0 <= x[1:2] <= 2.0, start = -1.5)
        @variable(model, -1.0 <= u[1:2] <= 1.0)
        @constraint(model, ∂(x[1]) == u[1])
        @constraint(model, ∂(x[2]) == u[2])
        return @constraint(model, x in Final(ball))
    end
    @test problem.target_set === ball
end

@testset "unconstrained final coordinates fall back to the variable bounds" begin
    # FIXED-L7 (Phase 3). Before Phase 3 the unconstrained coordinate contributed ±Inf, from
    # which `UT.box` built a NaN radius and threw inside LazySets.
    problem = lowered_problem() do model
        @variable(model, -1.0 <= x[1:2] <= 1.0, start = 0.0)
        @variable(model, -1.0 <= u[1:2] <= 1.0)
        @constraint(model, ∂(x[1]) == u[1])
        @constraint(model, ∂(x[2]) == u[2])
        # Only x[1] gets a target.
        return @constraint(model, final(x[1]) in MOI.Interval(-0.5, 0.5))
    end

    @test LazySets.low(problem.target_set, 1) == -0.5
    @test LazySets.high(problem.target_set, 1) == 0.5
    @test LazySets.low(problem.target_set, 2) == -1.0
    @test LazySets.high(problem.target_set, 2) == 1.0
end

@testset "horizon: seconds in continuous time, steps in discrete time" begin
    target = UT.box(SVector(-0.5), SVector(0.5))

    # No horizon ⇒ infinite.
    infinite = lowered_problem() do model
        integrator!(model)
        return @constraint(model, state(model) in Final(target))
    end
    @test infinite.time === PR.Infinity()

    # Reach is a "within at most T" spec, so 1.0 s at a 0.3 s step rounds *down*.
    reach = lowered_problem(; time_step = 0.3) do model
        integrator!(model)
        @constraint(model, state(model) in Final(target))
        return set_attribute(model, "horizon", 1.0)
    end
    @test reach.time == 3

    # Safety is "for at least T", so the same horizon rounds *up*.
    safety = lowered_problem(; time_step = 0.3) do model
        integrator!(model)
        @constraint(model, state(model) in Always(UT.box(SVector(-1.8), SVector(1.8))))
        return set_attribute(model, "horizon", 1.0)
    end
    @test safety.time == 4

    # A discrete-time model counts steps directly, with no `time_step` involved.
    discrete = lowered_problem() do model
        @variable(model, -2.0 <= x <= 2.0, start = -1.5)
        @variable(model, -1.0 <= u <= 1.0)
        @constraint(model, Δ(x) == x + u)
        @constraint(model, state(model) in Final(target))
        return set_attribute(model, "horizon", 7.0)
    end
    @test discrete.time == 7

    # A continuous-time horizon without a time step cannot be converted.
    @test_throws ErrorException lowered_problem() do model
        integrator!(model)
        @constraint(model, state(model) in Final(target))
        return set_attribute(model, "horizon", 1.0)
    end
end

@testset "validation" begin
    target = UT.box(SVector(-0.5), SVector(0.5))

    # FIXED-L9 (Phase 3): an objective is rejected instead of being silently dropped.
    err = try
        lowered_problem() do model
            integrator!(model)
            @constraint(model, state(model) in Final(target))
            return @objective(model, Min, model[:u])
        end
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("transition_cost", err.msg)

    # An unbounded variable cannot be discretized; the error names it rather than surfacing
    # as a NaN assertion from LazySets.
    err = try
        lowered_problem() do model
            @variable(model, x >= -2.0, start = -1.5)      # no upper bound
            @variable(model, -1.0 <= u <= 1.0)
            @constraint(model, ∂(x) == u)
            return @constraint(model, state(model) in Final(target))
        end
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("x", err.msg)
    @test occursin("Unbounded", err.msg)

    # A specification set must span the whole state vector.
    @test_throws ErrorException lowered_problem() do model
        @variable(model, -1.0 <= x[1:2] <= 1.0, start = 0.0)
        @variable(model, -1.0 <= u[1:2] <= 1.0)
        @constraint(model, ∂(x[1]) == u[1])
        @constraint(model, ∂(x[2]) == u[2])
        return @constraint(model, [x[1]] in Final(UT.box(SVector(-0.5), SVector(0.5))))
    end

    # Two sets of the same kind are ambiguous.
    @test_throws ErrorException lowered_problem() do model
        integrator!(model)
        @constraint(model, state(model) in Final(target))
        return @constraint(model, state(model) in Final(UT.box(SVector(0.0), SVector(1.0))))
    end
end

@testset "end-to-end: a safety problem solved from JuMP" begin
    # Safety was unreachable from the JuMP front-end before Phase 3.
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -2.0 <= x <= 2.0, start = 0.0)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == u)
    @constraint(model, state(model) in Always(UT.box(SVector(-1.5), SVector(1.5))))

    set_attribute(model, "time_step", 0.5)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "print_level", 0)

    optimize!(model)

    @test get_attribute(model, "concrete_problem") isa PR.SafetyProblem
    @test is_solved_and_feasible(model)

    traj = Dionysos.simulate(model, SVector(0.0); nsteps = 20)
    safe_set = get_attribute(model, "concrete_problem").safe_set
    @test all(x -> x ∈ safe_set, ST.states(traj))
end

@testset "end-to-end: a reach-avoid model avoids the unsafe region" begin
    # `ẋ = u` on [-2, 2], driving to a target around the origin. Whether the `Always` set is
    # present must change the *synthesized controller*, not just the lowered problem: with it,
    # states outside the safe set are no longer controllable, even though the target is still
    # perfectly reachable from them.
    function reach_avoid(safe)
        model = direct_model(Dionysos.Optimizer())
        @variable(model, -2.0 <= x <= 2.0, start = -1.5)
        @variable(model, -1.0 <= u <= 1.0)
        @constraint(model, ∂(x) == u)
        @constraint(model, state(model) in Final(UT.box(SVector(-0.25), SVector(0.25))))
        safe === nothing || @constraint(model, state(model) in Always(safe))

        set_attribute(model, "time_step", 0.5)
        set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
        set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
        set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
        set_attribute(model, "print_level", 0)

        optimize!(model)
        return model, get_attribute(model, "concrete_controller")
    end

    open_model, open_ctrl = reach_avoid(nothing)
    @test is_solved_and_feasible(open_model)
    @test controller_admissible(open_ctrl, SVector(-1.5))
    @test controller_admissible(open_ctrl, SVector(1.5))   # reachable from the right too

    fenced, fenced_ctrl = reach_avoid(UT.box(SVector(-1.8), SVector(0.5)))
    @test get_attribute(fenced, "concrete_problem").safe_set !== nothing
    @test is_solved_and_feasible(fenced)                   # the start is still safe
    @test controller_admissible(fenced_ctrl, SVector(-1.5))
    # x = 1.5 lies outside the safe set, so it is no longer part of the winning region.
    @test !controller_admissible(fenced_ctrl, SVector(1.5))
end

@testset "end-to-end: an abstraction-only model" begin
    # No specification: build the abstraction and hand it back.
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0)
    @variable(model, -1.0 <= u <= 1.0)
    @constraint(model, ∂(x) == u)

    set_attribute(model, "time_step", 0.5)
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "print_level", 0)

    optimize!(model)

    abstract_system = get_attribute(model, "abstract_system")
    @test abstract_system !== nothing
    @test SY.get_n_state(abstract_system) > 0
    # Nothing was synthesized, so there is no controller — and `simulate` says so plainly.
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.NO_SOLUTION
    @test_throws ErrorException Dionysos.simulate(model, SVector(0.0))
end

end # module TestMain
