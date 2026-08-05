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

# Phase 6: clocks. A clock is not new syntax — it is a state whose dynamics is the constant 1
# (running) or 0 (frozen). Its presence turns the augmented state from `(x, mode)` into
# `(x, t, mode)` and lets a target carry a time window.

const WR = Dionysos.Wrapper

# The thermostat again, now with a clock: the timer runs while the heater is off and freezes
# while it is on, so the model can ask to reach a temperature within a time window.
function clocked_model(; Ta = 18.0, α = 0.1, β = 2.0, tmax = 2.0)
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)

    @variable(model, 17.0 <= T <= 21.0, start = 18.0)
    @variable(model, 0.0 <= u <= 1.0)
    @variable(model, 0.0 <= t <= tmax)

    off = @mode(model, off)
    on = @mode(model, on)

    @constraint(off, ∂(T) == -α * (T - Ta))
    @constraint(on, ∂(T) == -α * (T - Ta) + β * u)
    @constraint(off, ∂(t) == 1)     # the timer runs while the heater is off…
    @constraint(on, ∂(t) == 0)      # …and is frozen while it is on

    @constraint(off, u == 0.0)
    @constraint(on, 0.5 <= u <= 1.0)

    add_transition!(model, off => on) do tr
        return @constraint(tr, T <= 19.0)
    end
    add_transition!(model, on => off) do tr
        return @constraint(tr, T >= 20.0)
    end

    return model, off, on
end

function lowered(model)
    opt = backend(model)
    WR.infer_roles!(opt.ir)
    return WR.build_problem(opt.ir, opt.dynamics_backend), opt
end

@testset "a constant-rate state is recognised as the clock" begin
    model, off, on = clocked_model()
    @constraint(off, [model[:T]] in Final(UT.box(SVector(19.5), SVector(21.0))))

    opt = backend(model)
    WR.infer_roles!(opt.ir)

    clock = WR.clock_index(opt.ir)
    @test clock !== nothing
    @test opt.ir.variables[clock].name == "t"
    # The clock is a state of the hybrid system, but not a coordinate of the physical `x`.
    @test WR.state_indices(opt.ir) == [1]
    @test WR.input_indices(opt.ir) == [2]
end

@testset "a clock makes each mode a physical system paired with its time axis" begin
    model, off, on = clocked_model()
    @constraint(off, [model[:T]] in Final(UT.box(SVector(19.5), SVector(21.0))))

    problem, _ = lowered(model)
    hs = problem.system

    for k in 1:2
        @test HybridSystems.mode(hs, k) isa ST.VectorContinuousSystem
    end

    # The timer runs in mode 1 (A = 1) and is frozen in mode 2 (A = 0).
    @test HybridSystems.mode(hs, 1).systems[2].A[1, 1] == 1.0
    @test HybridSystems.mode(hs, 2).systems[2].A[1, 1] == 0.0
    # …over the clock's declared domain.
    @test LazySets.high(MS.stateset(HybridSystems.mode(hs, 1).systems[2]), 1) == 2.0
end

@testset "guards and resets live in the augmented (x, t) space" begin
    model, off, on = clocked_model()
    @constraint(off, [model[:T]] in Final(UT.box(SVector(19.5), SVector(21.0))))

    problem, _ = lowered(model)
    hs = problem.system

    guard = MS.stateset(
        HybridSystems.resetmap(hs, first(HybridSystems.transitions(hs.automaton))),
    )
    # Two dimensions now: the temperature bound from the guard, and the whole time domain.
    @test LazySets.dim(guard) == 2
    @test LazySets.high(guard, 1) == 19.0
    @test LazySets.low(guard, 2) == 0.0
    @test LazySets.high(guard, 2) == 2.0
end

@testset "a clock reset restarts the timer" begin
    model, off, on = clocked_model()
    @constraint(off, [model[:T]] in Final(UT.box(SVector(19.5), SVector(21.0))))

    # A second on → off transition, distinct from the one the helper declared: several
    # transitions may share a (source, target) pair, and each keeps its own guard and reset.
    add_transition!(model, on => off) do tr
        @constraint(tr, model[:T] >= 20.5)
        return @constraint(tr, Δ(model[:t]) == 0)     # restart the timer on this switch
    end

    problem, _ = lowered(model)
    hs = problem.system
    resets =
        [HybridSystems.resetmap(hs, tr) for tr in HybridSystems.transitions(hs.automaton)]
    # The transition carrying `Δ(t) == 0` maps (T, t) to (T, 0); the others are the identity.
    outputs = [MS.apply(r, [20.6, 1.7]) for r in resets]
    @test any(o -> o ≈ [20.6, 0.0], outputs)
    @test any(o -> o ≈ [20.6, 1.7], outputs)
end

@testset "a target on the clock becomes a time window" begin
    model, off, on = clocked_model()
    target = UT.box(SVector(19.5), SVector(21.0))
    @constraint(off, [model[:T]] in Final(target))
    @constraint(off, final(model[:t]) in MOI.Interval(0.5, 1.5))

    problem, _ = lowered(model)

    @test problem.target_set isa PR.HybridSpec
    spec = problem.target_set.per_mode[1]
    @test spec isa PR.TimedSpec
    @test spec.tmin == 0.5
    @test spec.tmax == 1.5

    # Membership is over the augmented (x, t, mode) state.
    @test PR.satisfies(problem.target_set, SVector(20.0), 1.0, 1)
    @test !PR.satisfies(problem.target_set, SVector(20.0), 1.9, 1)   # outside the window
    @test !PR.satisfies(problem.target_set, SVector(18.0), 1.0, 1)   # outside the target

    # The initial state gains a time coordinate.
    @test problem.initial_set == (SVector(18.0), 0.0, 1)
end

@testset "end-to-end: a clocked hybrid model solved from JuMP" begin
    model, off, on = clocked_model(; tmax = 1.0)
    target = UT.box(SVector(19.0), SVector(21.0))
    @constraint(on, [model[:T]] in Final(target))

    for m in (off, on)
        set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
        set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
        set_attribute(m, "time_step", 0.5)
        set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
        set_attribute(m, "jacobian_bound", u -> SMatrix{1, 1}(-0.1))
    end

    optimize!(model)

    abstract_system = get_attribute(model, "abstract_system")
    @test abstract_system isa SY.HybridSymbolicModel
    # Every mode is clock-lifted, so the abstract states carry a time coordinate.
    q = first(SY.enum_states(abstract_system))
    @test length(SY.get_concrete_state(abstract_system, q)) == 3
    @test get_attribute(model, "concrete_controller") !== nothing
end

end # module TestMain
