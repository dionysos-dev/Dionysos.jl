module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Spot

# The accepting set of the co-safe reduction is soundness-critical under a folded (universal)
# run: an over-large set silently verifies violating states, where synthesis would at least fail
# visibly. These tests pin the finite-trace acceptance — a state accepts iff the run stopped
# there, padded with the empty valuation forever, is an accepted word — that replaced the
# absorbing-states heuristic, on the formulas that separate the two.

@testset "finite-trace acceptance is computed, not guessed" begin
    # The flagship's formula over the labels its quotient actually emits — observations are
    # mutually exclusive, so conjunctions never occur. Only the monitor's true-sink has every
    # obligation discharged.
    φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
    S = Dionysos.spot_stepper(φ)
    labels = [(), (:D,), (:R1,), (:R2,), (:R3,)]
    @test OPDS.accepting_states(S, labels) == Set([2])

    # Acceptance is per-state under the stop-and-pad reading, so the alphabet only prunes
    # reachability: the free-powerset judgement agrees with the emitted-labels one here.
    @test OPDS.accepting_states(S) == Set([2])

    # The reach-avoid shape is the reason the reading is finite-trace: no finite prefix
    # certifies the infinite `G`, but the "goal seen, hazard never" state accepts once the run
    # stops there — and the initial state does not.
    S_ra = Dionysos.spot_stepper(ltl"F(goal) & G(!hazard)")
    acc_ra = OPDS.accepting_states(S_ra, [(), (:goal,), (:hazard,)])
    @test !isempty(acc_ra)
    @test S_ra.qa0 ∉ acc_ra

    # The degenerate end of that spectrum: G(!a) alone is satisfied by the empty continuation
    # from its initial state — pure safety, no finite-trace content. Declaring everything
    # accepting would be unsound under verification; the computation must refuse instead.
    S_safety = Dionysos.spot_stepper(ltl"G(!a)")
    @test_throws ErrorException OPDS.accepting_states(S_safety, [(), (:a,)])

    # A monitor that states its accepting set explicitly is taken at its word, with or without
    # alphabet information.
    monitor = OPDS.FunctionMonitor(1, Set([2]), (qa, ap) -> qa)
    @test OPDS.accepting_states(monitor) == Set([2])
    @test OPDS.accepting_states(monitor, [()]) == Set([2])
end

end # module TestMain
