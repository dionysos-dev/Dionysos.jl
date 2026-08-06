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

# Phase 5: hybrid models written from JuMP — modes, transitions, guards, reset maps — checked
# against the hand-built `HybridSystem` the solver has always consumed.

const WR = Dionysos.Wrapper

# The thermostat: heater off / heater on, switching on temperature thresholds.
function thermostat_model(; Ta = 18.0, α = 0.1, β = 2.0)
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)

    @variable(model, 17.0 <= T <= 25.0, start = 18.0)
    @variable(model, 0.0 <= u <= 1.0)

    @mode(model, off)
    @mode(model, on)

    @constraint(off, ∂(T) == -α * (T - Ta))
    @constraint(on, ∂(T) == -α * (T - Ta) + β * u)

    @constraint(off, u == 0.0)              # per-mode input set
    @constraint(on, 0.2 <= u <= 1.0)

    add_transition!(model, off => on) do t
        return @constraint(t, T <= 19.0)
    end
    add_transition!(model, on => off) do t
        return @constraint(t, T >= 21.0)
    end

    return model, off, on
end

function lowered(model)
    opt = backend(model)
    return WR.lower(opt), opt
end

@testset "modes and transitions lower to a HybridSystem" begin
    model, off, on = thermostat_model()
    @constraint(off, [model[:T]] in Final(UT.box(SVector(21.0), SVector(23.0))))
    @constraint(on, [model[:T]] in Final(UT.box(SVector(21.0), SVector(23.0))))

    problem, _ = lowered(model)
    hs = problem.system

    @test hs isa HybridSystems.HybridSystem
    @test HybridSystems.nmodes(hs.automaton) == 2
    @test HybridSystems.ntransitions(hs.automaton) == 2

    # Each mode is a plain physical system: the augmented state is (x, mode), no clock.
    for k in 1:2
        @test HybridSystems.mode(hs, k) isa MS.ConstrainedBlackBoxControlContinuousSystem
    end

    # Per-mode input sets: the heater is off in mode 1 and throttled in mode 2. The bounds go
    # through JuMP's affine normalisation, so they are only exact to floating point.
    @test LazySets.high(MS.inputset(HybridSystems.mode(hs, 1)), 1) ≈ 0.0 atol = 1e-12
    @test LazySets.low(MS.inputset(HybridSystems.mode(hs, 2)), 1) ≈ 0.2

    # Per-mode dynamics: at T = Ta the off mode is at rest, the on mode is heating.
    @test MS.mapping(HybridSystems.mode(hs, 1))(SVector(18.0), SVector(0.0))[1] ≈ 0.0
    @test MS.mapping(HybridSystems.mode(hs, 2))(SVector(18.0), SVector(1.0))[1] ≈ 2.0
end

@testset "guards and reset maps" begin
    model, off, on = thermostat_model()
    @constraint(off, [model[:T]] in Final(UT.box(SVector(21.0), SVector(23.0))))

    problem, _ = lowered(model)
    hs = problem.system

    transitions = collect(HybridSystems.transitions(hs.automaton))
    by_source = Dict(HybridSystems.source(hs.automaton, t) => t for t in transitions)

    # off → on is enabled on T ≤ 19, intersected with the mode's own state set.
    guard_off = MS.stateset(HybridSystems.resetmap(hs, by_source[1]))
    @test LazySets.low(guard_off, 1) == 17.0
    @test LazySets.high(guard_off, 1) == 19.0

    # on → off is enabled on T ≥ 21.
    guard_on = MS.stateset(HybridSystems.resetmap(hs, by_source[2]))
    @test LazySets.low(guard_on, 1) == 21.0
    @test LazySets.high(guard_on, 1) == 25.0

    # With no `Δ` written on the transition the reset is the identity.
    @test MS.apply(HybridSystems.resetmap(hs, by_source[1]), SVector(18.5)) == SVector(18.5)
end

@testset "an explicit reset map is compiled" begin
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    @variable(model, -2.0 <= x <= 2.0, start = -1.0)
    @variable(model, -1.0 <= u <= 1.0)

    @mode(model, a)
    @mode(model, b)
    @constraint(a, ∂(x) == u)
    @constraint(b, ∂(x) == u)
    @constraint(b, [x] in Final(UT.box(SVector(-0.2), SVector(0.2))))

    add_transition!(model, a => b) do t
        @constraint(t, x >= 1.0)
        return @constraint(t, Δ(x) == -x)     # bounce off the boundary
    end

    problem, _ = lowered(model)
    hs = problem.system
    reset = HybridSystems.resetmap(hs, first(HybridSystems.transitions(hs.automaton)))
    @test MS.apply(reset, SVector(1.5))[1] ≈ -1.5
end

@testset "modes, transitions and scoped constraints display readably" begin
    # A scope subtypes `JuMP.AbstractModel`, so JuMP's generic `show` would try to summarise it
    # as a model and ask for `objective_sense`. Displaying one must work: it happens whenever a
    # user types `off` at the REPL, and whenever Documenter renders an example block — a docs
    # build failed on exactly this.
    #
    # `Model(…)` rather than `direct_model`, because that is what the documentation and the
    # README use, and its caching layer is what lets a constraint be printed back.
    model = Model(Dionysos.Optimizer)
    @variable(model, 17.0 <= T <= 25.0)
    @variable(model, 0.0 <= u <= 1.0)
    @mode(model, off)
    @mode(model, on)
    @constraint(off, ∂(T) == -0.1 * (T - 18.0))
    @constraint(on, ∂(T) == -0.1 * (T - 18.0) + 2.0 * u)

    transition = add_transition!(model, on => off) do t
        return @constraint(t, T >= 22.0)
    end

    @test sprint(show, MIME"text/plain"(), off) == "Mode(:off)"
    @test occursin("mode 2 => mode 1", sprint(show, MIME"text/plain"(), transition))

    # `print` needs its own method: JuMP defines `Base.print(::IO, ::AbstractModel)`, which is
    # more specific than the generic `Base.print` that falls back to `show`. Without it,
    # `println(off)` and `"$off"` reach JuMP's model summary and throw.
    @test sprint(print, off) == "Mode(:off)"
    @test "$off" == "Mode(:off)"
    @test occursin("mode 2 => mode 1", "$transition")

    # A scoped constraint shows its scope, not the wrapper's internal set types.
    scoped = @constraint(off, T <= 24.0)
    text = sprint(show, MIME"text/plain"(), scoped)
    @test occursin("[mode 1]", text)
    @test !occursin("ScopedSet", text)

    spec = @constraint(on, [T] in Final(UT.box(SVector(21.0), SVector(23.0))))
    spec_text = sprint(show, MIME"text/plain"(), spec)
    @test occursin("Final(", spec_text)
    @test !occursin("ScopedVectorSet", spec_text)
end

@testset "a hybrid reach-avoid model keeps its Always set" begin
    # `Final` and `Always` on the same hybrid model. The hybrid path has no state space to fold
    # an `Always` set into, so before `OptimalControlProblem` gained a `safe_set` it was simply
    # dropped here and the avoid part of the specification went unsynthesized.
    function thermostat_reach_avoid(; safe_low = nothing)
        model, off, on = thermostat_model()          # starts at T = 18
        target = UT.box(SVector(21.0), SVector(23.0))
        for m in (off, on)
            @constraint(m, [model[:T]] in Final(target))
            safe_low === nothing || @constraint(
                m,
                [model[:T]] in Always(UT.box(SVector(safe_low), SVector(25.0)))
            )
            set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
            set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
            set_attribute(m, "time_step", 0.5)
            set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
            set_attribute(m, "jacobian_bound", u -> SMatrix{1, 1}(-0.1))
        end
        return model
    end

    open_problem, _ = lowered(thermostat_reach_avoid())
    @test open_problem.safe_set === nothing

    problem = WR.lower(backend(thermostat_reach_avoid(; safe_low = 20.0)))
    @test problem isa PR.OptimalControlProblem
    @test problem.safe_set isa PR.HybridSpec
    @test PR.satisfies(problem.safe_set, SVector(22.0), 1)
    @test !PR.satisfies(problem.safe_set, SVector(18.0), 1)   # below the safe band

    # The model starts at T = 18, so a safe band that excludes it has no solution at all —
    # which is itself the proof that the set is enforced rather than dropped. Without the
    # `Always` constraint the very same model solves.
    fenced = thermostat_reach_avoid(; safe_low = 20.0)
    optimize!(fenced)
    @test !is_solved_and_feasible(fenced)

    open_model = thermostat_reach_avoid()
    optimize!(open_model)
    @test is_solved_and_feasible(open_model)

    # Widen the band to cover the start and it is feasible again, with the safe set reaching
    # the abstract problem and genuinely restricting it.
    reachable = thermostat_reach_avoid(; safe_low = 17.5)
    optimize!(reachable)
    @test is_solved_and_feasible(reachable)

    abstract_problem = backend(reachable).inner.control_solver.abstract_problem
    @test abstract_problem.safe_set !== nothing
    @test !isempty(abstract_problem.safe_set)
    nstates = SY.get_n_state(SY.get_automaton(get_attribute(reachable, "abstract_system")))
    @test length(abstract_problem.safe_set) < nstates
end

@testset "end-to-end: a half-space guard survives the abstraction" begin
    # A guard that is not a box reaches `get_states_from_set(source_model, guard, INNER)`. A
    # lazy `Intersection` enumerates there; an `IntersectionArray` would not, because the
    # discretisation resolves that one by computing a concrete `intersection`, which has no
    # method for most set pairs. `general_sets.jl` checks the lowered guard; this checks that
    # the abstraction can actually be built from it.
    model = direct_model(Dionysos.Optimizer())
    set_attribute(model, "print_level", 0)
    @variable(model, -2.0 <= x <= 2.0, start = -1.5)
    @variable(model, -2.0 <= y <= 2.0, start = -1.5)
    @variable(model, -1.0 <= u[1:2] <= 1.0)

    @mode(model, a)
    @mode(model, b)
    for m in (a, b)
        @constraint(m, ∂(x) == u[1])
        @constraint(m, ∂(y) == u[2])
        set_attribute(m, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))
        set_attribute(m, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0)))
        set_attribute(m, "time_step", 0.5)
        set_attribute(m, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    end
    @constraint(b, [x, y] in Final(UT.box(SVector(1.0, 1.0), SVector(1.5, 1.5))))

    add_transition!(model, a => b) do t
        return @constraint(t, x + y >= 0.0)      # a half-space, not a box
    end

    optimize!(model)

    @test get_attribute(model, "abstract_system") isa SY.HybridSymbolicModel
    # Switch transitions were built from the non-box guard, so the modes are connected.
    automaton = SY.get_automaton(get_attribute(model, "abstract_system"))
    @test SY.get_n_state(automaton) > 0
    @test !isempty(collect(SY.enum_transitions(automaton)))
end

@testset "@mode binds the name in the calling scope, like @variable" begin
    # No assignment needed: `@mode(model, cruise)` is enough, and `cruise = @mode(...)` is
    # merely redundant. Pinned because the documentation and every example rely on it.
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0)
    @variable(model, -1.0 <= u <= 1.0)

    @mode(model, cruise)
    @test cruise isa WR.Mode
    @test cruise.name === :cruise
    @test model[:cruise] === cruise      # registered in the object dictionary too

    @mode(model, coast)
    @test coast.id == cruise.id + 1      # ids keep advancing

    # The macro still returns the mode, so the assigned form keeps working.
    braking = @mode(model, braking)
    @test braking === model[:braking]
end

@testset "a transition without a guard is rejected" begin
    model = direct_model(Dionysos.Optimizer())
    @variable(model, -1.0 <= x <= 1.0)
    @variable(model, -1.0 <= u <= 1.0)
    @mode(model, a)
    @mode(model, b)
    @constraint(a, ∂(x) == u)
    @constraint(b, ∂(x) == u)

    @test_throws ErrorException add_transition!(model, a => b) do t
        return nothing
    end
end

@testset "the hybrid solver is selected, and modes carry their own options" begin
    model, off, on = thermostat_model()
    target = UT.box(SVector(21.0), SVector(23.0))
    @constraint(off, [model[:T]] in Final(target))
    @constraint(on, [model[:T]] in Final(target))

    grid = MP.GridFree(SVector(0.0), SVector(0.5))
    for m in (off, on)
        set_attribute(m, "state_grid", grid)
        set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
        set_attribute(m, "time_step", 0.5)
        set_attribute(m, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
    end

    opt = backend(model)
    problem = WR.lower(opt)
    @test WR.is_hybrid(opt.ir)
    @test WR.mode_ids(opt.ir) == [1, 2]
    # The options landed on the modes, not on the model-level solver.
    @test any(p -> first(p) == "state_grid", opt.ir.modes[1].attributes)

    @test WR.select_solver(problem.system, problem) === AB.HybridSystemAbstraction.Optimizer

    # The specification became a mode-indexed spec over the augmented state.
    @test problem.target_set isa PR.HybridSpec
    @test PR.satisfies(problem.target_set, SVector(22.0), 2)
    @test !PR.satisfies(problem.target_set, SVector(24.0), 2)
    # The initial state is the augmented (x, mode) pair the hybrid solver expects.
    @test problem.initial_set == (SVector(18.0), 1)
end

@testset "end-to-end: a hybrid model solved from JuMP" begin
    model, off, on = thermostat_model()
    target = UT.box(SVector(20.5), SVector(23.0))
    @constraint(off, [model[:T]] in Final(target))
    @constraint(on, [model[:T]] in Final(target))

    for m in (off, on)
        set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
        set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
        set_attribute(m, "time_step", 0.5)
        set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
        set_attribute(m, "jacobian_bound", u -> SMatrix{1, 1}(-0.1))
    end

    optimize!(model)

    @test backend(model).inner isa AB.HybridSystemAbstraction.Optimizer
    @test get_attribute(model, "abstract_system") isa SY.HybridSymbolicModel
    @test get_attribute(model, "concrete_controller") !== nothing

    traj = Dionysos.simulate(model, (SVector(18.0), 1); nsteps = 60)
    @test !isempty(ST.states(traj))
    @test first(ST.modes(traj)) == 1
end

end # module TestMain
