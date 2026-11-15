using Test

if !isdefined(Main, :LTL2Automata)
    include(joinpath(@__DIR__, "..", "..", "src", "LTL", "LTL2Automata.jl"))
end
const TestL2A = Main.LTL2Automata

@testset "infinite-horizon Büchi synthesis" begin
    # G(yellow): the sole accepting state has a transition only on yellow.
    always_yellow = TestL2A.Automaton(
        ["yellow"], 0, 1, true,
        Dict{Tuple{Int,UInt64},Int}((0, UInt64(1)) => 0),
        BitSet([0]),
    )

    safe_loop = TestL2A.TSFromAdj([[[1]]], [1])
    W, policy, next_q, _ = TestL2A.synthesize_buchi_controller(
        safe_loop, always_yellow, [Set([:yellow])]; semantics=:robust
    )
    @test W == BitVector([true])
    @test policy == [1]
    @test !isempty(TestL2A.post(safe_loop, 1, policy[1]))

    deadlock = TestL2A.TSFromAdj([[Int[]]], [1])
    Wdead, policy_dead, _, _ = TestL2A.synthesize_buchi_controller(
        deadlock, always_yellow, [Set([:yellow])]; semantics=:robust
    )
    @test Wdead == BitVector([false])
    @test policy_dead == [0]

    unsafe_loop = TestL2A.TSFromAdj([[[1]]], [1])
    Wunsafe, _, next_unsafe, _ = TestL2A.synthesize_buchi_controller(
        unsafe_loop, always_yellow, [Set{Symbol}()]; semantics=:robust
    )
    @test Wunsafe == BitVector([false])
    @test next_unsafe[1, 1] == TestL2A.I_INVALID_Q
    _, _, _, unsafe_status = TestL2A.simulate_buchi_policy(
        unsafe_loop, always_yellow, [Set{Symbol}()], [0], next_unsafe;
        x0=1, q0=always_yellow.init, consume_initial_label=true,
        stop_at_accept=false,
    )
    @test unsafe_status == :invalid_initial_transition

    # GF(yellow): q=0 is accepting; q=1 records a non-yellow suffix.
    infinitely_often_yellow = TestL2A.Automaton(
        ["yellow"], 0, 2, true,
        Dict{Tuple{Int,UInt64},Int}(
            (0, UInt64(1)) => 0,
            (0, UInt64(0)) => 1,
            (1, UInt64(1)) => 0,
            (1, UInt64(0)) => 1,
        ),
        BitSet([0]),
    )

    alternating = TestL2A.TSFromAdj([[[2]], [[1]]], [1])
    alternating_labels = [Set{Symbol}(), Set([:yellow])]
    Walt, palt, _, _ = TestL2A.synthesize_buchi_controller(
        alternating, infinitely_often_yellow, alternating_labels; semantics=:robust
    )
    @test all(Walt)
    @test all(>(0), palt)

    never_yellow = TestL2A.TSFromAdj([[[1]]], [1])
    Wnever, pnever, _, _ = TestL2A.synthesize_buchi_controller(
        never_yellow, infinitely_often_yellow, [Set{Symbol}()]; semantics=:robust
    )
    @test !any(Wnever)
    @test all(==(0), pnever)
end


@testset "threaded Büchi predecessors across packed-word boundaries" begin
    # The product sizes deliberately leave thread partitions unaligned with
    # BitVector's 64-bit words. Every plant state has its own self-loop, so
    # G(yellow) wins exactly at yellow-labelled states, without a solver oracle.
    always_yellow = TestL2A.Automaton(
        ["yellow"], 0, 1, true,
        Dict{Tuple{Int,UInt64},Int}((0, UInt64(1)) => 0), BitSet([0]),
    )
    for n in (257, 1025, 4101)
        ts = TestL2A.TSFromAdj([[[x], Int[]] for x in 1:n], collect(1:n))
        for pattern in (:all_winning, :alternating, :mixed)
            expected = BitVector([
                pattern == :all_winning ||
                (pattern == :alternating ? isodd(x) : mod(x, 7) != 0)
                for x in 1:n
            ])
            labels = [expected[x] ? Set([:yellow]) : Set{Symbol}() for x in 1:n]
            for repetition in 1:16
                W, policy, _, _ = TestL2A.synthesize_buchi_controller(
                    ts, always_yellow, labels; semantics=:robust, max_fp_iters=20,
                )
                @test W == expected
                @test policy == Int.(expected)
            end
        end
    end
end
