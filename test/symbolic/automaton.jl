module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import HybridSystems

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

# Sample transitions in the `add_transition!(autom, source, target, symbol)`
# order. Chosen to exercise nondeterminism on (source = 5, symbol = 6), a
# self-loop (5 → 5), and several distinct (source, symbol) branches.
const SAMPLE = [(5, 9, 6), (5, 8, 6), (5, 3, 7), (8, 3, 6), (5, 5, 6), (8, 3, 7)]

# The same transitions in the stored `(target, source, symbol)` order.
const SAMPLE_STORED =
    Set([(9, 5, 6), (8, 5, 6), (3, 5, 7), (3, 8, 6), (5, 5, 6), (3, 8, 7)])

const NSTATES = 10
const NSYMBOLS = 11

function _populate!(autom)
    for (source, target, symbol) in SAMPLE
        SY.add_transition!(autom, source, target, symbol)
    end
    return autom
end

# A subtype that overrides nothing, to check the abstract interface stubs error.
struct _UnimplementedAutom <: SY.AbstractAutomatonList end

# ---------------------------------------------------------------------------
# Interface conformance — runs identically for every concrete implementation.
# All membership checks compare sets, so successor/predecessor *ordering* is
# deliberately left out of the contract.
# ---------------------------------------------------------------------------

function run_automaton_interface_tests(make::Function)
    label = nameof(typeof(make(NSTATES, NSYMBOLS)))
    @testset "$label" begin
        @testset "construction & sizes" begin
            autom = make(NSTATES, NSYMBOLS)
            @test SY.get_n_state(autom) == NSTATES
            @test SY.get_n_input(autom) == NSYMBOLS
            @test SY.enum_states(autom) == 1:NSTATES
            @test SY.enum_inputs(autom) == 1:NSYMBOLS
            @test HybridSystems.ntransitions(autom) == 0
        end

        @testset "empty-automaton queries" begin
            autom = make(NSTATES, NSYMBOLS)
            @test isempty(SY.post(autom, 5, 6))
            @test isempty(collect(SY.pre(autom, 3)))
            @test SY.is_deterministic(autom)
            @test isempty(SY.nondeterminism_counts(autom))
            @test SY.count_self_loops(autom) == 0
        end

        @testset "add_transition! / enum / post / pre" begin
            autom = _populate!(make(NSTATES, NSYMBOLS))
            @test HybridSystems.ntransitions(autom) == length(SAMPLE)
            @test Set(SY.enum_transitions(autom)) == SAMPLE_STORED
            @test Set(SY.post(autom, 5, 6)) == Set([9, 8, 5])
            @test Set(collect(SY.pre(autom, 3))) == Set([(5, 7), (8, 6), (8, 7)])
            @test isempty(SY.post(autom, 5, 1))       # unknown (source, symbol)
            @test isempty(collect(SY.pre(autom, 4)))   # unknown target
        end

        @testset "compute_post! appends in place" begin
            autom = _populate!(make(NSTATES, NSYMBOLS))
            out = Int[7]                 # pre-seeded: compute_post! must append, not clobber
            SY.compute_post!(out, autom, 5, 6)
            @test Set(out) == Set([7, 9, 8, 5])
            SY.compute_post!(out, autom, 8, 7)
            @test Set(out) == Set([7, 9, 8, 5, 3])
        end

        @testset "bulk add_transitions! matches individual" begin
            bulk = make(NSTATES, NSYMBOLS)
            SY.add_transitions!(bulk, collect(SAMPLE_STORED))
            @test HybridSystems.ntransitions(bulk) == length(SAMPLE)
            @test Set(SY.enum_transitions(bulk)) == SAMPLE_STORED
            @test Set(SY.post(bulk, 5, 6)) == Set([9, 8, 5])
        end

        @testset "determinism analytics" begin
            det = make(NSTATES, NSYMBOLS)
            SY.add_transition!(det, 5, 9, 6)
            @test SY.is_deterministic(det)

            autom = _populate!(make(NSTATES, NSYMBOLS))
            @test !SY.is_deterministic(autom)
            counts = SY.nondeterminism_counts(autom)
            @test sum(counts) == length(SAMPLE)
            @test maximum(counts) == 3            # (source = 5, symbol = 6) has 3 successors
            @test SY.count_self_loops(autom) == 1  # the 5 → 5 transition
        end

        @testset "add_state!" begin
            autom = make(NSTATES, NSYMBOLS)
            new = HybridSystems.add_state!(autom)
            @test new == NSTATES + 1
            @test SY.get_n_state(autom) == NSTATES + 1
            @test SY.enum_states(autom) == 1:(NSTATES + 1)

            SY.add_transition!(autom, new, 1, 1)   # new → 1
            SY.add_transition!(autom, 1, new, 1)   # 1 → new
            @test Set(SY.post(autom, new, 1)) == Set([1])
            @test Set(collect(SY.pre(autom, new))) == Set([(1, 1)])
            @test Set(SY.post(autom, 1, 1)) == Set([new])
            @test Set(collect(SY.pre(autom, 1))) == Set([(new, 1)])
        end

        @testset "finalize! preserves the transition relation" begin
            autom = _populate!(make(NSTATES, NSYMBOLS))
            SY.add_transition!(autom, 5, 9, 6)     # duplicate of an existing transition
            @test SY.finalize!(autom) === autom
            @test Set(SY.post(autom, 5, 6)) == Set([9, 8, 5])
            @test Set(collect(SY.pre(autom, 3))) == Set([(5, 7), (8, 6), (8, 7)])
        end

        @testset "empty! clears transitions, keeps states" begin
            autom = _populate!(make(NSTATES, NSYMBOLS))
            Base.empty!(autom)
            @test HybridSystems.ntransitions(autom) == 0
            @test isempty(SY.post(autom, 5, 6))
            @test isempty(collect(SY.pre(autom, 3)))
            @test SY.count_self_loops(autom) == 0
            @test SY.get_n_state(autom) == NSTATES  # empty! drops transitions, not states
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Run the shared suite against every implementation.
# ---------------------------------------------------------------------------

run_automaton_interface_tests((n, m) -> SY.SortedAutomatonList(n, m))
run_automaton_interface_tests((n, m) -> SY.IndexedAutomatonList(n, m))
run_automaton_interface_tests((n, m) -> SY.FastIndexedAutomatonList(n, m))

# ---------------------------------------------------------------------------
# Implementation-specific behaviour.
# ---------------------------------------------------------------------------

@testset "FastIndexedAutomatonList: finalize! deduplicates" begin
    autom = _populate!(SY.FastIndexedAutomatonList(NSTATES, NSYMBOLS))
    SY.add_transition!(autom, 5, 9, 6)  # duplicate of (target = 9, source = 5, symbol = 6)
    @test HybridSystems.ntransitions(autom) == length(SAMPLE) + 1
    @test SY.post(autom, 5, 6) == [9, 8, 5, 9]  # duplicate present before finalize!
    SY.finalize!(autom)
    @test HybridSystems.ntransitions(autom) == length(SAMPLE)
    @test SY.post(autom, 5, 6) == [9, 8, 5]     # deduplicated in place
end

@testset "AbstractAutomatonList interface stubs error" begin
    autom = _UnimplementedAutom()
    @test_throws ErrorException SY.get_n_state(autom)
    @test_throws ErrorException SY.get_n_input(autom)
    @test_throws ErrorException SY.enum_transitions(autom)
    @test_throws ErrorException SY.add_transition!(autom, 1, 2, 3)
    @test_throws ErrorException SY.pre(autom, 1)
    @test_throws ErrorException SY.post(autom, 1, 2)
    @test_throws ErrorException Base.empty!(autom)
end

# ---------------------------------------------------------------------------
# Closed-loop simulation directly on an automaton (`closed_loop.jl`). Only the
# `post` map is implementation-specific, and it is exercised above for every
# backend, so one representative automaton suffices here.
# ---------------------------------------------------------------------------

# A static controller that commands input `1` on every state in `1:nstates`.
function _constant_input_controller(nstates::Int)
    table = ST.ControlTable(nstates)
    for q in 1:nstates
        ST.add_control!(table, q, 1)
    end
    return ST.DiscreteStaticController(Set(1:nstates), table, false)
end

@testset "get_closed_loop_trajectory: deterministic path stops at a dead end" begin
    # Path 1 → 2 → 3 → 4 under the single input symbol; state 4 has no successor.
    autom = SY.IndexedAutomatonList(4, 1)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 2, 3, 1)
    SY.add_transition!(autom, 3, 4, 1)

    traj = SY.get_closed_loop_trajectory(autom, _constant_input_controller(4), 1, 10)
    # `first` successor each step; halts when `post(4, 1)` is empty.
    @test collect(ST.states(traj)) == [1, 2, 3, 4]
end

@testset "get_closed_loop_trajectory: randomize_post picks an admissible successor" begin
    # (source = 1, symbol = 1) is nondeterministic: successors {1, 2}.
    autom = SY.IndexedAutomatonList(2, 1)
    SY.add_transition!(autom, 1, 1, 1)  # self-loop
    SY.add_transition!(autom, 1, 2, 1)

    traj = SY.get_closed_loop_trajectory(
        autom,
        _constant_input_controller(2),
        1,
        5;
        randomize_post = true,
    )
    visited = collect(ST.states(traj))
    @test visited[1] == 1
    @test all(q -> q in (1, 2), visited)  # every chosen successor is admissible
end

@testset "get_closed_loop_trajectory: stopping predicate halts early" begin
    autom = SY.IndexedAutomatonList(4, 1)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 2, 3, 1)
    SY.add_transition!(autom, 3, 4, 1)

    traj = SY.get_closed_loop_trajectory(
        autom,
        _constant_input_controller(4),
        1,
        10;
        stopping = q -> q == 3,
    )
    @test collect(ST.states(traj)) == [1, 2, 3]
end

end # module
