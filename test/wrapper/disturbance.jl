module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Dionysos
using JuMP
using Symbolics
using MathOptSymbolicAD
import MathOptInterface as MOI
import LazySets

const WR = Dionysos.Wrapper

function lowered_problem(build!)
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    build!(model)
    opt = backend(model)
    return WR.lower(opt), opt
end

@testset "a declared disturbance lowers to a Noisy system" begin
    problem, _ = lowered_problem() do model
        @variable(model, -1.0 <= x <= 1.0, start = -0.75)
        @variable(model, -1.0 <= u <= 1.0)
        @variable(model, -0.1 <= w <= 0.1)
        set_role!(w, Dionysos.DISTURBANCE)
        @constraint(model, ∂(x) == -x + u + w)
        return @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    end

    system = problem.system
    @test system isa MathematicalSystems.NoisyConstrainedBlackBoxControlContinuousSystem
    # The disturbance set is the variable's bounds — the bounds *are* the declaration.
    W = MathematicalSystems.noiseset(system)
    @test LazySets.low(W) == [-0.1] && LazySets.high(W) == [0.1]
    # The compiled dynamics take the three-place signature and mean what the model wrote.
    @test system.f([0.5], [0.2], [0.1]) ≈ [-0.5 + 0.2 + 0.1]

    # The same model with `w` undeclared reads it as an input — the roles are what separate
    # robust synthesis from synthesis over a bigger input space.
    plain, _ = lowered_problem() do model
        @variable(model, -1.0 <= x <= 1.0, start = -0.75)
        @variable(model, -1.0 <= u <= 1.0)
        @variable(model, -0.1 <= w <= 0.1)
        @constraint(model, ∂(x) == -x + u + w)
        return @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    end
    @test plain.system isa MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem
    @test MathematicalSystems.inputdim(plain.system) == 2

    # Discrete time follows the same rule.
    discrete, _ = lowered_problem() do model
        @variable(model, -1.0 <= x <= 1.0, start = -0.75)
        @variable(model, -1.0 <= u <= 1.0)
        @variable(model, -0.1 <= w <= 0.1)
        set_role!(w, Dionysos.DISTURBANCE)
        @constraint(model, Δ(x) == 0.5 * x + u + w)
        return @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    end
    @test discrete.system isa
          MathematicalSystems.NoisyConstrainedBlackBoxControlDiscreteSystem
end

@testset "what a disturbance must not be" begin
    # A driven variable is a state, whoever perturbs it.
    @test_throws ErrorException lowered_problem() do model
        @variable(model, -1.0 <= x <= 1.0, start = 0.0)
        @variable(model, -0.1 <= w <= 0.1)
        set_role!(w, Dionysos.DISTURBANCE)
        @constraint(model, ∂(w) == -w)
        @constraint(model, ∂(x) == -x + w)
        return @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    end

    # The hybrid lowering does not thread a disturbance yet, and must say so.
    @test_throws ErrorException lowered_problem() do model
        @variable(model, -1.0 <= x <= 1.0, start = 0.0)
        @variable(model, -0.1 <= w <= 0.1)
        set_role!(w, Dionysos.DISTURBANCE)
        @mode(model, m1)
        @constraint(model, m1, ∂(x) == -x + w)
        return @constraint(model, final(x) in MOI.Interval(-0.5, 0.5))
    end
end

@testset "the symbolic extension derives the noisy bounds" begin
    X = LazySets.Hyperrectangle(; low = [-2.0], high = [2.0])
    U = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])

    # Non-additive shape (noisedim ≠ statedim): the additive reading is unavailable, so the
    # bound zᵢ ≥ sup |fᵢ(x,u,w) − fᵢ(x,u,w_c)| is derived by tracing the difference — the
    # dynamics cancel symbolically and the 2w₁ term bounds exactly.
    skewed = MathematicalSystems.NoisyConstrainedBlackBoxControlContinuousSystem(
        (x, u, w) -> SVector(-x[1] + u[1] + 2.0 * w[1] + 0.0 * w[2]),
        1,
        1,
        2,
        X,
        U,
        LazySets.Hyperrectangle(; low = [-0.1, -0.1], high = [0.1, 0.1]),
    )
    @test ST.compute_noise_bound(skewed) ≈ [0.2]

    # With Symbolics loaded the growth-bound constructor no longer needs the bounds spelled
    # out: the Jacobian is bounded with w ranged over W, and z is derived.
    approx = ST.ContinuousTimeGrowthBound(skewed)
    @test approx isa ST.ContinuousTimeGrowthBound

    # The Jacobian bound of a noisy system ranges the disturbance too: ∂f/∂x = -1 + w here,
    # and freezing w at its centre would report 1.0 where the sound bound over W is 1.5.
    wobbly = MathematicalSystems.NoisyConstrainedBlackBoxControlContinuousSystem(
        (x, u, w) -> SVector((-1.0 + w[1]) * x[1] + u[1]),
        1,
        1,
        1,
        X,
        U,
        LazySets.Hyperrectangle(; low = [-0.5], high = [0.5]),
    )
    L = ST.compute_jacobian_bound(wobbly)
    @test L(SVector(0.0))[1, 1] ≈ -0.5
end

@testset "end to end: the front-end synthesizes robustly" begin
    model = Model(Dionysos.Optimizer)
    @variable(model, -2.0 <= x <= 2.0, start = -1.0)
    @variable(model, -1.0 <= u <= 1.0)
    @variable(model, -0.2 <= w <= 0.2)
    set_role!(w, Dionysos.DISTURBANCE)
    @constraint(model, ∂(x) == -x + u + w)
    @constraint(model, final(x) in MOI.Interval(-0.3, 0.3))

    set_attribute(model, "time_step", 0.5)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
    set_attribute(model, "print_level", 0)
    optimize!(model)

    @test get_attribute(model, "concrete_controller") !== nothing
    concrete_system = get_attribute(model, "concrete_problem").system
    @test concrete_system isa
          MathematicalSystems.NoisyConstrainedBlackBoxControlContinuousSystem
end

end # module TestMain
