module TestBoundedInputVariation

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI

# ----------------------------------------------------------------------------
# 3×3 grid walk (states 1..9 row-major), moves UP/DOWN/LEFT/RIGHT/STAY.
# Input "vectors" give the slew-rate distance: with Δ = 1 turns and
# start/stop are allowed but direct reversals (UP↔DOWN, LEFT↔RIGHT) are not.
# ----------------------------------------------------------------------------

const UP, DOWN, LEFT, RIGHT, STAY = 1, 2, 3, 4, 5
const MOVE = Dict(
    UP => SVector(0.0, 1.0),
    DOWN => SVector(0.0, -1.0),
    LEFT => SVector(-1.0, 0.0),
    RIGHT => SVector(1.0, 0.0),
    STAY => SVector(0.0, 0.0),
)
input_distance(u1::Int, u2::Int) = maximum(abs.(MOVE[u1] - MOVE[u2]))

function grid_automaton()
    autom = SY.IndexedAutomatonList(9, 5)
    for q in 1:9
        r, c = divrem(q - 1, 3) .+ (1, 1)
        r > 1 && SY.add_transition!(autom, q, q - 3, UP)
        r < 3 && SY.add_transition!(autom, q, q + 3, DOWN)
        c > 1 && SY.add_transition!(autom, q, q - 1, LEFT)
        c < 3 && SY.add_transition!(autom, q, q + 1, RIGHT)
        SY.add_transition!(autom, q, q, STAY)
    end
    return autom
end

# Deterministic single successor.
step(autom, q, u) = first(SY.post(autom, q, u))

# Run the dynamic controller on the automaton, returning the visited states and
# played inputs (protocol order matches the closed-loop engine: output then
# update, both on the same memory and measurement).
function simulate(autom, ctrl, q0; nstep = 20, stop = q -> false)
    qs, us = [q0], Int[]
    mem = ST.initial_state(ctrl)
    q = q0
    for _ in 1:nstep
        stop(q) && break
        u = ST.output_control(ctrl, mem, q)
        u === nothing && break
        mem = ST.update_state(ctrl, mem, q)
        q = step(autom, q, u)
        push!(us, u)
        push!(qs, q)
    end
    return qs, us
end

@testset "grid walk with turn-limited inputs" begin
    autom = grid_automaton()

    constraint = OPDS.BoundedInputVariation(input_distance, 1.0)
    ctrl, controllable, _, value_tab = OPDS.compute_bounded_input_variation_controller(
        autom,
        [9],
        constraint;
        initial_set = [1],
    )

    # Two RIGHTs and two DOWNs in some interleaved order: turns cost nothing
    # extra under Δ = 1, so the constrained optimum equals the unconstrained 4.
    @test value_tab[1] == 4.0
    @test 1 in controllable

    qs, us = simulate(autom, ctrl, 1; stop = q -> q == 9)
    @test qs[end] == 9
    @test length(us) == 4
    @test all(input_distance(us[k], us[k + 1]) <= 1.0 for k in 1:(length(us) - 1))
end

@testset "Δ too small: no turn possible" begin
    autom = grid_automaton()

    constraint = OPDS.BoundedInputVariation(input_distance, 0.5)
    _, controllable, _, value_tab = OPDS.compute_bounded_input_variation_controller(
        autom,
        [9],
        constraint;
        initial_set = [1],
    )
    # 1 → 9 needs a turn (RIGHT then DOWN), forbidden at Δ = 0.5.
    @test value_tab[1] == Inf
    @test 1 ∉ controllable

    # Straight-line target on the same row stays feasible.
    _, controllable3, _, value_tab3 = OPDS.compute_bounded_input_variation_controller(
        autom,
        [3],
        constraint;
        initial_set = [1],
    )
    @test value_tab3[1] == 2.0
    @test 1 in controllable3
end

@testset "boundary inputs: ramp down into the target" begin
    autom = grid_automaton()

    # The final input entering the target must be within Δ of STAY, which every
    # unit move satisfies at Δ = 1 — but not at Δ = 0.9.
    constraint = OPDS.BoundedInputVariation(input_distance, 0.9; target_input = STAY)
    _, controllable, _, _ = OPDS.compute_bounded_input_variation_controller(
        autom,
        [9],
        constraint;
        initial_set = [1],
    )
    @test 1 ∉ controllable
end

# ----------------------------------------------------------------------------
# The pitfall that motivates pair-keyed values: a per-state controller table
# commits state 2 to its locally cheapest input c2, which is incompatible with
# arriving via a — the executed pair then violates the constraint. The
# pair-keyed controller plays c1 when arriving via a and c2 when starting
# fresh.
# ----------------------------------------------------------------------------
@testset "memory-dependent input choice (per-state table pitfall)" begin
    a, c1, c2 = 1, 2, 3
    level = Dict(a => 0.0, c1 => 0.0, c2 => 10.0)
    d(u1, u2) = abs(level[u1] - level[u2])

    autom = SY.IndexedAutomatonList(3, 3)
    SY.add_transition!(autom, 1, 2, a)
    SY.add_transition!(autom, 2, 3, c1)
    SY.add_transition!(autom, 2, 3, c2)
    cost(q, u) = (q, u) == (2, c1) ? 10.0 : (q, u) == (2, c2) ? 2.0 : 1.0

    constraint = OPDS.BoundedInputVariation(d, 1.0)
    ctrl, controllable, _, value_tab = OPDS.compute_bounded_input_variation_controller(
        autom,
        [3],
        constraint;
        initial_set = [1],
        cost_function = cost,
    )

    @test 1 in controllable
    @test value_tab[1] == 1.0 + 10.0 # forced through the expensive c1
    @test value_tab[2] == 2.0        # fresh start at 2 may use the cheap c2

    # Fresh at 2: cheapest input. Arriving via a: the only compatible one.
    @test ST.output_control(ctrl, 0, 2) == c2
    @test ST.output_control(ctrl, a, 2) == c1

    qs, us = simulate(autom, ctrl, 1; stop = q -> q == 3)
    @test qs == [1, 2, 3]
    @test us == [a, c1]
end

@testset "nondeterministic automata are rejected" begin
    autom = SY.IndexedAutomatonList(3, 1)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 1, 3, 1) # two successors for (1, 1)
    constraint = OPDS.BoundedInputVariation((u1, u2) -> 0.0, 1.0)
    @test_throws ErrorException OPDS.compute_bounded_input_variation_controller(
        autom,
        [3],
        constraint,
    )
end

@testset "optimizer wiring" begin
    autom = grid_automaton()
    problem = PR.OptimalControlProblem(autom, [1], [9], nothing, nothing, PR.Infinity())

    optimizer = OPDS.OptimizerOptimalControlProblem()
    optimizer.problem = problem
    optimizer.print_level = 0
    optimizer.bounded_input_variation = OPDS.BoundedInputVariation(input_distance, 1.0)
    MOI.optimize!(optimizer)

    @test optimizer.success
    @test optimizer.controller isa ST.DiscreteDynamicController
    @test optimizer.value_function(1) == 4.0

    qs, _ = simulate(autom, optimizer.controller, 1; stop = q -> q == 9)
    @test qs[end] == 9
end

end # module
