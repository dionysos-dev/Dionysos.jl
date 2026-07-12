module TestClockLift

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathematicalSystems as MS

# Small hand-built base abstraction: 2 states, 1 input, 2 transitions (1→2, 2→2).
function build_base()
    Xmap = MP.ExplicitGridMapping(MP.GridFree(SVector(0.0), SVector(1.0)))
    MP.add_pos!(Xmap, (0,))   # q1 at x = 0.0
    MP.add_pos!(Xmap, (1,))   # q2 at x = 1.0
    Umap = MP.ExplicitGridMapping(MP.GridFree(SVector(0.0), SVector(1.0)))
    MP.add_pos!(Umap, (0,))   # u1
    base = SY.SymbolicModelList(Xmap, Umap)
    SY.add_transitions!(base.autom, [(2, 1, 1), (2, 2, 1)])  # (q′, q, u)
    return base
end

active_clock(tmax, tstep) = SY.ClockAbstraction(
    MS.ConstrainedLinearContinuousSystem([1.0;;], UT.box([0.0], [tmax])),
    tstep,
)
frozen_clock(tmax) = SY.ClockAbstraction(
    MS.ConstrainedLinearContinuousSystem([0.0;;], UT.box([0.0], [tmax])),
    1.0,
)

@testset "ClockLift - active clock (product structure)" begin
    base = build_base()
    clock = active_clock(2.0, 1.0)          # tsteps = [0, 1, 2]  (3 steps)
    @test clock.is_active
    @test length(clock.tsteps) == 3

    m = SY.lift(SY.ClockLift(clock), base)

    # State count = base states × time steps; inputs pass through.
    @test SY.get_n_state(m) == 2 * 3
    @test SY.get_n_input(m) == SY.get_n_input(base)
    @test SY.get_state_dim(m) == SY.get_state_dim(base) + 1

    # Each base transition is replicated across (ntime - 1) advancing steps.
    @test SY.get_n_transitions(m) == 2 * 2

    # Every lifted transition advances the clock by exactly one step (p → p+1).
    for tr in SY.enum_transitions(m)
        src = SY.get_concrete_state(m, SY.transition_source(tr))
        tgt = SY.get_concrete_state(m, SY.transition_target(tr))
        @test tgt[end] ≈ src[end] + 1.0
    end

    # Concrete ↔ abstract round-trip over the augmented coordinate [x; t].
    for id in SY.enum_states(m)
        @test SY.get_abstract_state(m, SY.get_concrete_state(m, id)) == id
    end
end

@testset "ClockLift - get_states_from_set projects the [x; t] box" begin
    base = build_base()
    m = SY.lift(SY.ClockLift(active_clock(2.0, 1.0)), base)

    # Box covers both cells (x ∈ {0,1}) and time steps p=1,2 (t ∈ {0,1}); excludes t=2.
    box = UT.box(SVector(-0.5, -0.5), SVector(1.5, 1.5))
    states = SY.get_states_from_set(m, box, MP.INNER)

    @test length(states) == 2 * 2
    # Soundness: every returned state's concrete coordinate lies in the query box.
    for id in states
        @test SY.get_concrete_state(m, id) ∈ box
    end
end

@testset "ClockLift - frozen clock keeps a single time slice" begin
    base = build_base()
    clock = frozen_clock(2.0)
    @test !clock.is_active

    m = SY.lift(SY.ClockLift(clock), base)

    @test SY.get_n_state(m) == SY.get_n_state(base)   # no time factor
    @test SY.get_n_transitions(m) == SY.get_n_transitions(base)
    for id in SY.enum_states(m)
        @test SY.get_concrete_state(m, id)[end] == 0.0       # frozen at t = 0
        @test SY.get_abstract_state(m, SY.get_concrete_state(m, id)) == id
    end
end

end # module
