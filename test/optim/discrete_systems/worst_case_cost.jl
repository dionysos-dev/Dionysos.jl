module TestMain

import Dionysos
import Random
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

# `compute_optimal_controller`: the min-max Dijkstra search with a general stage cost.

# ----------------------------------------------------------------------------
# A cheaper path found after a dearer one must replace it in the controller. From 1, `d` reaches
# the target 3 at cost 5; `c` goes through 2 at cost 1 + 1, but 2 is popped after 1 is relaxed.
# ----------------------------------------------------------------------------
@testset "a strict improvement replaces the control (deterministic)" begin
    c, d = 1, 2
    autom = SY.IndexedAutomatonList(3, 2)
    SY.add_transition!(autom, 1, 3, d)
    SY.add_transition!(autom, 1, 2, c)
    SY.add_transition!(autom, 2, 3, c)
    cost(q, u) = (q, u) == (1, d) ? 5.0 : 1.0
    ctrl, _, _, value = OPDS.compute_optimal_controller(autom, [3]; cost_function = cost)
    @test value[1] == 2.0
    @test ST.output_control(ctrl, nothing, 1) == c
end

# ----------------------------------------------------------------------------
# Non-deterministic: `a` from 1 lands on 2 (value 1) or 3 (value 4), so it costs 1 + 4 = 5, which
# beats the direct `b` at cost 10. `b` completes first, at the pop of the target, so `a` must
# replace it, with the worst successor taken from the pop order.
# ----------------------------------------------------------------------------
@testset "min-max with a non-uniform cost" begin
    a, b = 1, 2
    autom = SY.IndexedAutomatonList(4, 2)
    SY.add_transition!(autom, 1, 2, a)
    SY.add_transition!(autom, 1, 3, a)
    SY.add_transition!(autom, 1, 4, b)
    SY.add_transition!(autom, 2, 4, a)
    SY.add_transition!(autom, 3, 4, a)
    cost(q, u) = (q, u) == (1, b) ? 10.0 : (q, u) == (3, a) ? 4.0 : 1.0
    ctrl, _, _, value = OPDS.compute_optimal_controller(autom, [4]; cost_function = cost)
    @test value[2] == 1.0
    @test value[3] == 4.0
    @test value[1] == 5.0
    @test ST.output_control(ctrl, nothing, 1) == a
end

# ----------------------------------------------------------------------------
# On random non-deterministic automata with a unit cost, the search agrees with the breadth-first
# solver on every value.
# ----------------------------------------------------------------------------
@testset "agrees with the uniform-cost solver" begin
    rng = Random.MersenneTwister(0)
    for _ in 1:20
        n, m = 30, 3
        autom = SY.IndexedAutomatonList(n, m)
        for q in 2:n, u in 1:m
            rand(rng) < 0.7 || continue
            for _ in 1:rand(rng, 1:3)
                SY.add_transition!(autom, q, rand(rng, 1:n), u)
            end
        end
        _, _, _, v_general =
            OPDS.compute_optimal_controller(autom, [1]; cost_function = (q, u) -> 1.0)
        _, _, _, v_uniform = OPDS.compute_worst_case_uniform_cost_controller(autom, [1])
        @test v_general == v_uniform
    end
end

end # module
