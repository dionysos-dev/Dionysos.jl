module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import HybridSystems

# The point of the view is that the *unmodified* solvers compute the robust game over it, so these
# tests call the same two entry points the synthesis tests call — `compute_largest_invariant_set`
# and `compute_worst_case_cost_controller` — and only ever vary the automaton handed to them.

@testset "construction and interface" begin
    autom = SY.IndexedAutomatonList(3, 2)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 1, 3, 2)
    SY.add_transition!(autom, 2, 3, 1)

    view = SY.FoldedAutomaton(autom)
    @test SY.get_n_state(view) == 3
    @test SY.get_n_input(view) == 1
    @test SY.inner(view) === autom
    @test SY.inner_symbols(view, 1) == [1, 2]

    # A faithful re-labelling: same triples with the symbol folded, duplicates preserved.
    @test sort(collect(SY.enum_transitions(view))) ==
          sort([(2, 1, 1), (3, 1, 1), (3, 2, 1)])
    @test sort(collect(SY.pre(view, 3))) == [(1, 1), (2, 1)]
    @test sort(SY.post(view, 1, 1)) == [2, 3]

    # `pre` and the counters' initialisation enumerate identically, so a fold that collapses two
    # parallel edges onto the same triple must yield it twice — deduplicating on one side only
    # would desynchronise the fixed-point counters.
    dup = SY.IndexedAutomatonList(2, 2)
    SY.add_transition!(dup, 1, 2, 1)
    SY.add_transition!(dup, 1, 2, 2)
    dup_view = SY.FoldedAutomaton(dup)
    @test collect(SY.pre(dup_view, 2)) == [(1, 1), (1, 1)]
    @test SY.post(dup_view, 1, 1) == [2, 2]

    # Read-only: mutation goes through the inner automaton, never the view.
    @test_throws ErrorException SY.add_transition!(view, 1, 2, 1)
    @test_throws ErrorException HybridSystems.add_state!(view)
    @test_throws ErrorException empty!(view)

    # Validation: every inner symbol mapped, every folded symbol inhabited.
    @test_throws ArgumentError SY.FoldedAutomaton(autom, [1])
    @test_throws ArgumentError SY.FoldedAutomaton(autom, [1, 3])
    @test_throws ArgumentError SY.FoldedAutomaton(autom, [0, 1])
end

@testset "safety: the folded view is the robust invariant" begin
    # Mode 1 keeps state 2 safe (2 → 1); mode 2 throws it away (2 → 3). A controller that picks
    # the mode saves 2; an environment that picks it kills 2.
    autom = SY.IndexedAutomatonList(3, 2)
    SY.add_transition!(autom, 1, 1, 1)
    SY.add_transition!(autom, 1, 1, 2)
    SY.add_transition!(autom, 2, 1, 1)
    SY.add_transition!(autom, 2, 3, 2)
    SY.add_transition!(autom, 3, 3, 1)
    SY.add_transition!(autom, 3, 3, 2)

    safelist = [1, 2]

    _, invariant_exists, _ = OPDS.compute_largest_invariant_set(autom, safelist)
    _, invariant_forall, _ =
        OPDS.compute_largest_invariant_set(SY.FoldedAutomaton(autom), safelist)

    @test sort(collect(invariant_exists)) == [1, 2]
    @test sort(collect(invariant_forall)) == [1]

    # One mode: nothing to fold, so bare and view must agree exactly.
    single = SY.IndexedAutomatonList(3, 1)
    SY.add_transition!(single, 1, 1, 1)
    SY.add_transition!(single, 2, 1, 1)
    SY.add_transition!(single, 3, 3, 1)
    _, inv_bare, _ = OPDS.compute_largest_invariant_set(single, safelist)
    _, inv_view, _ =
        OPDS.compute_largest_invariant_set(SY.FoldedAutomaton(single), safelist)
    @test sort(collect(inv_bare)) == sort(collect(inv_view))
end

@testset "reachability: synthesis ⊋ robust ⊋ verification on one automaton" begin
    # Four inner symbols are the pairs (u, w) with u the controller's half and w the
    # environment's: 1 = (u₁,w₁), 2 = (u₁,w₂), 3 = (u₂,w₁), 4 = (u₂,w₂).
    #
    # From state 1 every choice can be spoiled by the environment; from state 4 the choice u₂
    # cannot; state 2 is a sink off the target. So with target {3}:
    #     pure synthesis  (four symbols)  wins {1, 3, 4} — it also picks the w,
    #     robust          (fold by u)     wins {3, 4}    — u₂ survives both w at 4 only,
    #     verification    (fold all)      wins {3}       — nothing is chosen at all.
    autom = SY.IndexedAutomatonList(4, 4)
    SY.add_transition!(autom, 1, 3, 1)
    SY.add_transition!(autom, 1, 2, 2)
    SY.add_transition!(autom, 1, 3, 3)
    SY.add_transition!(autom, 1, 2, 4)
    for s in 1:4
        SY.add_transition!(autom, 2, 2, s)
        SY.add_transition!(autom, 3, 3, s)
    end
    SY.add_transition!(autom, 4, 2, 1)
    SY.add_transition!(autom, 4, 2, 2)
    SY.add_transition!(autom, 4, 3, 3)
    SY.add_transition!(autom, 4, 3, 4)

    targetlist = [3]

    # No `initial_set`: the default covers every state, so the fixed point runs to completion
    # instead of early-stopping the moment the initial states are won.
    solve(a) = OPDS.compute_worst_case_cost_controller(a, targetlist)

    _, win_exists, _, _ = solve(autom)
    contr_robust, win_robust, _, _ = solve(SY.FoldedAutomaton(autom, [1, 1, 2, 2]))
    _, win_forall, _, _ = solve(SY.FoldedAutomaton(autom))

    @test sort(collect(win_exists)) == [1, 3, 4]
    @test sort(collect(win_robust)) == [3, 4]
    @test sort(collect(win_forall)) == [3]

    # The controller speaks the folded alphabet: at state 4 it must have committed to u₂, the
    # choice whose every environment refinement reaches the target.
    @test only(contr_robust.controller_map(4)) == 2
end

@testset "complete_with_sink turns a missing environment move into a loss" begin
    # State 1 has a successor under mode 1 and *no transition* under mode 2. Folded as is, the
    # union quietly forgets that the environment can play mode 2 there, and the invariant keeps
    # state 1 — a verification the abstraction has no right to. The pessimistic completion routes
    # the missing move to a losing sink, and the claim disappears.
    autom = SY.IndexedAutomatonList(2, 2)
    SY.add_transition!(autom, 1, 1, 1)
    SY.add_transition!(autom, 2, 1, 1)
    SY.add_transition!(autom, 2, 1, 2)

    completed, sink, ncompleted = SY.complete_with_sink(autom)
    @test sink == 3
    @test ncompleted == 1
    @test SY.post(completed, 1, 2) == [sink]
    @test SY.post(completed, sink, 1) == [sink]
    @test SY.post(completed, sink, 2) == [sink]

    safelist = [1, 2]
    _, unsound, _ = OPDS.compute_largest_invariant_set(SY.FoldedAutomaton(autom), safelist)
    _, sound, _ =
        OPDS.compute_largest_invariant_set(SY.FoldedAutomaton(completed), safelist)

    @test 1 in unsound       # the vacuous win the completion exists to prevent
    @test isempty(sound)     # pessimistic: the unknown move is assumed to lose
end

end # module TestMain
