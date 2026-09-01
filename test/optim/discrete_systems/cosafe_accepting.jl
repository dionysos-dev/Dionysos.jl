module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Spot

# The accepting set of the co-safe reduction is soundness-critical under a folded (universal)
# run: an over-large set silently verifies violating states, where synthesis would at least fail
# visibly. These tests pin the exact good-prefix computation that replaced the absorbing-states
# heuristic, on the formulas that separate the two.

@testset "good-prefix acceptance is computed, not guessed" begin
    # The flagship's formula over the labels its quotient actually emits — observations are
    # mutually exclusive, so conjunctions never occur. The monitor's true-sink, and only it,
    # certifies a good prefix.
    φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
    S = Dionysos.spot_stepper(φ)
    labels = [(), (:D,), (:R1,), (:R2,), (:R3,)]
    @test OPDS.accepting_states(S, labels) == Set([2])

    # The heuristic's blind spot: G(!a) is a safety formula with no good prefix — no finite word
    # guarantees it. The old code read the missing `a`-edge as a self-loop and declared its
    # state accepting; the computation must refuse instead.
    S_safety = Dionysos.spot_stepper(ltl"G(!a)")
    @test_throws ErrorException OPDS.accepting_states(S_safety, [(), (:a,)])

    # The alphabet matters: judged against the free powerset, the flagship formula has no good
    # prefix (a state may lack edges for conjunctions like R2 ∧ D that the system cannot emit),
    # and refusing there is the conservative direction.
    @test_throws ErrorException OPDS.accepting_states(S)

    # A monitor that states its accepting set explicitly is taken at its word, with or without
    # alphabet information.
    monitor = OPDS.FunctionMonitor(1, Set([2]), (qa, ap) -> qa)
    @test OPDS.accepting_states(monitor) == Set([2])
    @test OPDS.accepting_states(monitor, [()]) == Set([2])
end

end # module TestMain
