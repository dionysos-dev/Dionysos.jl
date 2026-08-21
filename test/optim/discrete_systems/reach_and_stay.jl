module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI

# `◇□ T` on hand-built automata, where the winning set can be worked out by hand. The two
# readings of "and stay" differ, and the difference is the point of most of what follows.

function solve(
    autom,
    target,
    safe;
    initial = collect(1:SY.get_n_state(autom)),
    stay_on_first_entry = false,
    early_stop = false,
)
    optimizer = MOI.instantiate(OPDS.OptimizerReachAndStayProblem)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("problem"),
        PR.ReachAndStayProblem(
            autom,
            initial,
            target,
            safe,
            PR.Infinity();
            stay_on_first_entry = stay_on_first_entry,
        ),
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), early_stop)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)
    return (
        controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("controller")),
        winning = sort(
            collect(MOI.get(optimizer, MOI.RawOptimizerAttribute("winning_set"))),
        ),
    )
end

# 1 --u1--> 2 --u1--> 3 --u1--> 4 --u1--> 4
# The target is {1, 4}: 4 is invariant, 1 is in the target but its only move leaves it.
function departure_automaton()
    autom = SY.IndexedAutomatonList(4, 1)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 2, 3, 1)
    SY.add_transition!(autom, 3, 4, 1)
    SY.add_transition!(autom, 4, 4, 1)
    return autom
end

@testset "◇□ tolerates a finite departure from the target" begin
    # From 1 the run leaves the target through 2 and 3 and settles in 4 forever. That is
    # *eventually* always in the target, so 1 is winning.
    r = solve(departure_automaton(), [1, 4], collect(1:4))
    @test r.winning == [1, 2, 3, 4]

    # Every winning cell gets exactly one input, fixed when it is first won and never revised —
    # the memorylessness Li & Liu's algorithm is named for.
    for q in r.winning
        @test length(collect(r.controller.controller_map(q))) == 1
    end
end

@testset "stay-on-first-entry refuses that departure" begin
    # The same run enters the target at 1 and immediately leaves it, which this reading
    # forbids, so 1 drops out — and cannot be routed through either.
    r = solve(departure_automaton(), [1, 4], collect(1:4); stay_on_first_entry = true)
    @test r.winning == [2, 3, 4]
    @test !ST.is_defined(r.controller, nothing, 1)
end

@testset "the safe set blocks the approach" begin
    # State 3 is the only route from {1, 2} to the target, so removing it strands them.
    r = solve(departure_automaton(), [4], [1, 2, 4])
    @test r.winning == [4]

    r = solve(departure_automaton(), [4], collect(1:4))
    @test r.winning == [1, 2, 3, 4]
end

@testset "a target with no outgoing transition never stays" begin
    # 5 is in the target but the system has nowhere to go from it, so it is not invariant and
    # nothing may be routed through it either.
    autom = SY.IndexedAutomatonList(5, 1)
    SY.add_transition!(autom, 1, 5, 1)
    SY.add_transition!(autom, 2, 3, 1)
    SY.add_transition!(autom, 3, 3, 1)

    for first_entry in (false, true)
        r = solve(autom, [3, 5], collect(1:5); stay_on_first_entry = first_entry)
        @test !(5 in r.winning)
        @test !(1 in r.winning)      # its only move lands in 5
        @test 2 in r.winning && 3 in r.winning
    end
end

@testset "nondeterminism is resolved against the controller" begin
    # u1 from 1 may land in 2 (winning) or in 4 (a dead end), so it is not usable; u2 is.
    autom = SY.IndexedAutomatonList(4, 2)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 1, 4, 1)
    SY.add_transition!(autom, 1, 3, 2)
    SY.add_transition!(autom, 3, 3, 2)
    SY.add_transition!(autom, 2, 2, 1)

    r = solve(autom, [2, 3], collect(1:4))
    @test 1 in r.winning
    @test collect(r.controller.controller_map(1)) == [2]
    @test !(4 in r.winning)
end

@testset "a cell joining the core later is not relaxed twice" begin
    # The counters behind the attractor may be decremented once per transition, ever. State 2
    # is won on the way in during the first round and *then* joins the invariant core in the
    # second, which is the one shape that tempts an implementation to relax it again. Doing so
    # drops the count of `(1, u1)` to zero while its other successor, the dead end 3, is still
    # unwon — and state 1 is reported winning on an input that may strand the run.
    autom = SY.IndexedAutomatonList(4, 1)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 1, 3, 1)     # same input, so u1 may land in either
    SY.add_transition!(autom, 2, 4, 1)
    SY.add_transition!(autom, 4, 4, 1)
    # 3 has no transition at all

    r = solve(autom, [4], collect(1:4))
    @test r.winning == [2, 4]
    @test !(1 in r.winning)
    @test !(3 in r.winning)
end

@testset "early_stop returns a winning set that still covers the initial set" begin
    # A long chain: 1 → 2 → … → 12, with 12 invariant. Stopping as soon as state 8 is won must
    # not drop state 8, and must not claim more than the full solve does.
    n = 12
    autom = SY.IndexedAutomatonList(n, 1)
    for q in 1:(n - 1)
        SY.add_transition!(autom, q, q + 1, 1)
    end
    SY.add_transition!(autom, n, n, 1)

    full = solve(autom, [n], collect(1:n); initial = [8])
    early = solve(autom, [n], collect(1:n); initial = [8], early_stop = true)

    @test full.winning == collect(1:n)
    @test 8 in early.winning
    @test issubset(early.winning, full.winning)
end

end # module TestMain
