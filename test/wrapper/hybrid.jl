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
    @constraint(
        off,
        [model[:T]] in
        Final(LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0)))
    )
    @constraint(
        on,
        [model[:T]] in
        Final(LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0)))
    )

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
    @constraint(
        off,
        [model[:T]] in
        Final(LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0)))
    )

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
    @constraint(
        b,
        [x] in Final(LazySets.Hyperrectangle(; low = SVector(-0.2), high = SVector(0.2)))
    )

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

    spec = @constraint(
        on,
        [T] in Final(LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0)))
    )
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
        target = LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0))
        for m in (off, on)
            @constraint(m, [model[:T]] in Final(target))
            safe_low === nothing || @constraint(
                m,
                [model[:T]] in Always(
                    LazySets.Hyperrectangle(;
                        low = SVector(safe_low),
                        high = SVector(25.0),
                    ),
                )
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
    # Wide enough to contain whole cells: on a 0.5 grid an `INNER` target needs a cell centre at
    # least half a step inside it, and `[1.0, 1.5]²` — which this used to ask for — admits only a
    # centre at 1.25, which is not a grid point. The target discretized to nothing and the test
    # still passed, because it only checks that the abstraction was built.
    @constraint(
        b,
        [x, y] in Final(
            LazySets.Hyperrectangle(;
                low = SVector(0.75, 0.75),
                high = SVector(1.75, 1.75),
            ),
        )
    )

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
    target = LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0))
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
    # The initial condition is mode-indexed too, so it is a set over the augmented state rather
    # than the single `(x, mode)` point it used to collapse to.
    @test problem.initial_set isa PR.HybridSpec
    @test PR.satisfies(problem.initial_set, SVector(18.0), 1)
    @test !PR.satisfies(problem.initial_set, SVector(18.0), 2)
end

@testset "end-to-end: a hybrid model solved from JuMP" begin
    model, off, on = thermostat_model()
    target = LazySets.Hyperrectangle(; low = SVector(20.5), high = SVector(23.0))
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

# A switched system — the switch *is* the control, so no mode carries a continuous input. Its
# input space is zero-dimensional, which holds exactly one point: the action of leaving the
# switch alone. Nothing needs to be declared for it; the front-end supplies that grid.
@testset "a mode whose only control is the switch" begin
    safe = LazySets.Hyperrectangle(; low = SVector(-2.0), high = SVector(2.0))

    model = Model(Dionysos.Optimizer)
    @variable(model, -2.0 <= x <= 2.0)
    @mode(model, up)
    @mode(model, down)
    @constraint(up, ∂(x) == 1.0 - 0.5 * x)
    @constraint(down, ∂(x) == -1.0 - 0.5 * x)

    add_transition!(model, up => down) do t
        return @constraint(t, [x] in Guard(safe))
    end
    add_transition!(model, down => up) do t
        return @constraint(t, [x] in Guard(safe))
    end

    @constraint(up, [x] in Always(safe))
    @constraint(down, [x] in Always(safe))

    for m in (up, down)
        set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
        set_attribute(m, "time_step", 0.1)
        set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
        set_attribute(m, "jacobian_bound", u -> SMatrix{1, 1}(-0.5))
        set_attribute(m, "print_level", 0)
    end

    optimize!(model)                     # no `input_grid` set anywhere
    @test is_solved_and_feasible(model)

    # One continuous input per mode — the "hold the switch" action, without which the state
    # could not evolve — plus one switching input per transition.
    gim = get_attribute(model, "abstract_system").input_mapping
    @test gim.continuous_inputs == 2
    @test gim.switching_inputs == 2

    traj = Dionysos.simulate(model, (SVector(0.0), 1); nsteps = 40)
    @test all(x ∈ safe for x in ST.states(traj))
end

@testset "a mode without an `Always` set is unconstrained, not forbidden" begin
    # Writing `Always` on one mode used to leave every other mode out of the `HybridSpec`
    # entirely, and an absent mode does not satisfy the spec — so the specification silently
    # outlawed every mode it was not written on. A mode nobody constrained is now bounded by its
    # own state set instead.
    band = LazySets.Hyperrectangle(; low = SVector(18.0), high = SVector(24.0))

    model, off, _on = thermostat_model()
    @constraint(off, [model[:T]] in Always(band))
    problem, _ = lowered(model)

    @test problem isa PR.SafetyProblem
    @test sort(collect(keys(problem.safe_set.per_mode))) == [1, 2]
    @test PR.satisfies(problem.safe_set, SVector(20.0), 1)
    @test PR.satisfies(problem.safe_set, SVector(20.0), 2)     # was `false`
    # The mode that *was* constrained still is: the band, not its whole 17–25 range.
    @test !PR.satisfies(problem.safe_set, SVector(17.5), 1)
    # The unconstrained mode keeps its own bounds, so it is not unbounded either.
    @test !PR.satisfies(problem.safe_set, SVector(26.0), 2)
end

@testset "`EventuallyAlways` on a mode reaches the solver" begin
    # It used to be parsed onto the mode and then dropped: `build_hybrid_problem` read only
    # `FINAL` and `ALWAYS`, so a model carrying both `Always` and `EventuallyAlways` lowered to
    # the identical `SafetyProblem` and was solved for half its specification.
    band = LazySets.Hyperrectangle(; low = SVector(17.5), high = SVector(24.5))
    settle = LazySets.Hyperrectangle(; low = SVector(20.0), high = SVector(22.0))

    safety, _ = lowered(let (m, off, on) = thermostat_model()
        for k in (off, on)
            @constraint(k, [m[:T]] in Always(band))
        end
        m
    end)
    @test safety isa PR.SafetyProblem

    ras, _ = lowered(let (m, off, on) = thermostat_model()
        for k in (off, on)
            @constraint(k, [m[:T]] in Always(band))
            @constraint(k, [m[:T]] in EventuallyAlways(settle))
        end
        m
    end)
    @test ras isa PR.ReachAndStayProblem
    @test ras.target_set isa PR.HybridSpec
    @test ras.safe_set isa PR.HybridSpec
    @test PR.satisfies(ras.target_set, SVector(21.0), 1)
    @test !PR.satisfies(ras.target_set, SVector(23.0), 1)
    @test !ras.stay_on_first_entry

    # The flag travels from the marker to the problem.
    strict, _ = lowered(
        let (m, off, on) = thermostat_model()
            for k in (off, on)
                @constraint(k, [m[:T]] in EventuallyAlways(settle; stay_on_first_entry = true))
            end
            m
        end,
    )
    @test strict.stay_on_first_entry
    # With no `Always` the run is safe anywhere: the safe set is each mode's own state set.
    @test PR.satisfies(strict.safe_set, SVector(17.2), 1)

    # `◇□ T` already reaches `T`; asking for both is a contradiction worth naming.
    @test_throws ErrorException lowered(
        let (m, off, on) = thermostat_model()
            @constraint(off, [m[:T]] in Final(settle))
            @constraint(on, [m[:T]] in EventuallyAlways(settle))
            m
        end,
    )
end

@testset "a switch that discretizes to nothing fails, and losses are reported" begin
    # A guard no cell lies inside used to be a `@warn` and a `continue`: the modes ended up
    # disconnected, and synthesis then answered — often `OPTIMAL` — for a system the user had
    # not written. `INNER` inclusion needs a whole cell inside the guard, so an interval
    # narrower than one cell discretizes to nothing.
    function thin_guard(; width = 0.05)
        model = direct_model(Dionysos.Optimizer())
        set_attribute(model, "print_level", 0)
        @variable(model, 0.0 <= x <= 4.0, start = 0.5)
        @variable(model, -1.0 <= u <= 1.0)
        @mode(model, a)
        @mode(model, b)
        for m in (a, b)
            @constraint(m, ∂(x) == u)
            @constraint(
                m,
                [x] in
                Always(LazySets.Hyperrectangle(; low = SVector(0.0), high = SVector(4.0)))
            )
            set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
            set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
            set_attribute(m, "time_step", 0.5)
        end
        add_transition!(model, a => b) do t
            return @constraint(t, 2.0 <= x <= 2.0 + width)
        end
        return model
    end

    @test_throws ErrorException optimize!(thin_guard(; width = 0.05))

    # Wide enough to hold a cell, and the build report is clean.
    model = thin_guard(; width = 1.5)
    optimize!(model)
    report = SY.build_report(get_attribute(model, "abstract_system"))
    @test isempty(report)
    @test isempty(report.dropped_resets)
    @test isempty(report.inexact_resets)
end

@testset "a hybrid model starts from a set, not from its centre" begin
    # `Start(S)` used to collapse to `LazySets.center(S)`: feasibility was decided for one point
    # and every other state of the set the user declared went unsynthesized.
    function banded(; start_low = 17.5, start_high = 19.5)
        model, off, on = thermostat_model()
        @constraint(
            off,
            [model[:T]] in Start(
                LazySets.Hyperrectangle(;
                    low = SVector(start_low),
                    high = SVector(start_high),
                ),
            )
        )
        for m in (off, on)
            @constraint(
                m,
                [model[:T]] in Always(
                    LazySets.Hyperrectangle(; low = SVector(17.0), high = SVector(25.0)),
                )
            )
            set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
            set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
            set_attribute(m, "time_step", 0.5)
            set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
            set_attribute(m, "jacobian_bound", u -> SMatrix{1, 1}(-0.1))
        end
        return model
    end

    problem, _ = lowered(banded())
    @test problem.initial_set isa PR.HybridSpec
    @test collect(keys(problem.initial_set.per_mode)) == [1]   # the mode carrying `Start`

    model = banded()
    optimize!(model)
    @test is_solved_and_feasible(model)

    abstract_system = get_attribute(model, "abstract_system")
    q0s = SY.states_satisfying(abstract_system, WR.lower(backend(model)).initial_set)
    @test length(q0s) > 1                     # the whole 2 °C band, not its midpoint
    @test all(SY.get_concrete_state(abstract_system, q)[end] == 1 for q in q0s)

    # A start outside every mode's domain is named as such, rather than reported as an
    # infeasible control problem.
    @test_throws ErrorException optimize!(banded(; start_low = 40.0, start_high = 41.0))
end

@testset "end-to-end: hybrid reach-and-stay" begin
    # The band has to contain an equilibrium the input grid can actually hold. With `u` gridded
    # at 0.25 the heater's steady states are 18 + 20u = 23, 28, … °C, so 23 is the only one
    # inside the 17–25 domain: a band ending *at* 23 is one the abstraction cannot certify
    # staying in, and the winning set comes back empty.
    model, off, on = thermostat_model()
    settle = LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(25.0))
    for m in (off, on)
        @constraint(m, [model[:T]] in EventuallyAlways(settle))
        set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
        set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
        set_attribute(m, "time_step", 0.5)
        set_attribute(m, "approx_mode", AB.UniformGridAbstraction.GROWTH)
        set_attribute(m, "jacobian_bound", u -> SMatrix{1, 1}(-0.1))
    end

    optimize!(model)
    @test is_solved_and_feasible(model)

    winning = get_attribute(model, "winning_set")
    @test winning !== nothing
    @test !isempty(winning)

    # The run settles into the band and holds it. `stopping` is overridden because the default
    # stops on arrival, which is where a `◇□` run gets interesting.
    traj = Dionysos.simulate(model, (SVector(18.0), 1); nsteps = 120, stopping = _ -> false)
    states = collect(ST.states(traj))
    @test length(states) == 121          # the run is not cut short by an undefined control
    @test all(x ∈ settle for x in states[(end - 20):end])
end

end # module TestMain
