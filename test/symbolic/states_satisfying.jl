module TestStatesSatisfying

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathematicalSystems as MS
import LazySets

# Small hand-built base abstraction: 2 states, 1 input (as in the ClockLift test).
function build_base()
    Xmap = MP.ExplicitGridMapping(MP.GridFree(SVector(0.0), SVector(1.0)))
    MP.add_pos!(Xmap, (0,))   # q1 at x = 0.0
    MP.add_pos!(Xmap, (1,))   # q2 at x = 1.0
    Umap = MP.ExplicitGridMapping(MP.GridFree(SVector(0.0), SVector(1.0)))
    MP.add_pos!(Umap, (0,))
    base = SY.SymbolicModelList(Xmap, Umap)
    SY.add_transitions!(base.autom, [(2, 1, 1), (2, 2, 1)])
    return base
end

active_clock(tmax, tstep) = SY.ClockAbstraction(
    MS.ConstrainedLinearContinuousSystem(
        [1.0;;],
        LazySets.Hyperrectangle(; low = [0.0], high = [tmax]),
    ),
    tstep,
)

@testset "satisfies - concrete-point membership" begin
    box = LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))
    s = PR.StateSpec(box)
    @test PR.satisfies(s, SVector(0.0))
    @test !PR.satisfies(s, SVector(2.0))
    @test PR.satisfies(s, SVector(0.0), 5.0)          # plain spec ignores time

    ts = PR.TimedSpec(s, 1.0, 2.0)
    @test PR.satisfies(ts, SVector(0.0), 1.5)
    @test !PR.satisfies(ts, SVector(0.0), 0.5)        # outside time window
    @test !PR.satisfies(ts, SVector(2.0), 1.5)        # outside state set

    hs = PR.hybrid_reach_spec(
        [box],
        [LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0))],
        [2],
    )
    @test PR.satisfies(hs, SVector(0.0), 1.5, 2)
    @test !PR.satisfies(hs, SVector(0.0), 1.5, 1)     # mode not in spec
    @test !PR.satisfies(hs, SVector(0.0), 0.5, 2)     # outside time window
end

@testset "states_satisfying - clock-lifted model + lifted specs" begin
    m = SY.lift(SY.ClockLift(active_clock(2.0, 1.0)), build_base())  # tsteps [0,1,2]
    s = PR.StateSpec(LazySets.Hyperrectangle(; low = SVector(-0.5), high = SVector(1.5)))            # both cells

    # Base spec on a clock-lifted model: matching base states at every time index.
    all_states = SY.states_satisfying(m, s)
    @test length(all_states) == 2 * 3
    for id in all_states
        @test SVector(SY.get_concrete_state(m, id)[1]) ∈ s.set
    end

    # Timed spec: restrict to t ∈ [0,1] → time indices {1,2}.
    ts = PR.TimedSpec(s, 0.0, 1.0)
    timed_states = SY.states_satisfying(m, ts)
    @test length(timed_states) == 2 * 2
    for id in timed_states
        xt = SY.get_concrete_state(m, id)
        @test SVector(xt[1]) ∈ s.set
        @test 0.0 <= xt[end] <= 1.0
    end
end

end # module
