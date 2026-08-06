module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathematicalSystems as MS

# `GuardedResetMap` is the guard/reset pair stored on a `HybridSystems` transition. It replaces
# the per-problem `…ResetMap` structs that every hybrid model used to define for itself.

@testset "GuardedResetMap: guard is the state set, reset is the map" begin
    guard = UT.box(SVector(0.2), SVector(1.0))

    # The reset defaults to the identity — switching mode usually leaves the state alone.
    identity_map = ST.GuardedResetMap(guard)
    @test MS.stateset(identity_map) === guard
    @test MS.apply(identity_map, SVector(0.5)) == SVector(0.5)

    # A jump to a fixed point.
    fixed = ST.GuardedResetMap(guard, _ -> SVector(0.0))
    @test MS.apply(fixed, SVector(0.7)) == SVector(0.0)
    @test MS.stateset(fixed) === guard

    # A reset over the augmented `[x; t]` state, as a clock-lifted mode uses: reset x, clamp t.
    timed_guard = UT.box(SVector(0.0, 0.0), SVector(1.0, 5.0))
    timed = ST.GuardedResetMap(timed_guard, s -> vcat(0.0, max(2.0, s[end])))
    @test MS.apply(timed, [0.8, 1.0]) == [0.0, 2.0]
    @test MS.apply(timed, [0.8, 3.0]) == [0.0, 3.0]
end

@testset "GuardedResetMap drives a HybridSystem transition" begin
    import HybridSystems

    X = UT.box(SVector(-1.0), SVector(1.0))
    U = UT.box(SVector(-1.0), SVector(1.0))
    mode = MS.ConstrainedBlackBoxControlContinuousSystem((x, u) -> u, 1, 1, X, U)

    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    guard = UT.box(SVector(0.5), SVector(1.0))
    hs = HybridSystems.HybridSystem(
        automaton,
        [mode, mode],
        [ST.GuardedResetMap(guard)],
        [HybridSystems.AutonomousSwitching()],
    )

    transition = first(HybridSystems.transitions(hs.automaton))
    @test MS.stateset(HybridSystems.resetmap(hs, transition)) === guard
    @test MS.apply(HybridSystems.resetmap(hs, transition), SVector(0.6)) == SVector(0.6)
end

end # module TestMain
