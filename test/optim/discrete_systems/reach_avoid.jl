module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import MathOptInterface as MOI

# The `safe_set` of `PR.OptimalControlProblem`: the avoid part of a reach-avoid task, checked
# on hand-built automata so the expected controllable set can be worked out by hand.

# 1 --u1--> 2 --u1--> 3          the short way, through state 2
# 1 --u2--> 4 --u1--> 5 --u1--> 3   the long way round
function detour_automaton()
    autom = SY.IndexedAutomatonList(5, 2)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 2, 3, 1)
    SY.add_transition!(autom, 1, 4, 2)
    SY.add_transition!(autom, 4, 5, 1)
    SY.add_transition!(autom, 5, 3, 1)
    return autom
end

@testset "without a safe set the synthesis is unchanged" begin
    autom = detour_automaton()
    _, controllable, _, value =
        OPDS.compute_worst_case_cost_controller(autom, [3]; initial_set = [1])

    @test 1 in controllable
    @test 2 in controllable
    @test value[1] == 2                     # the short way, through state 2
end

@testset "an unsafe state is excluded, forcing the detour" begin
    autom = detour_automaton()
    controller, controllable, uncontrollable, value =
        OPDS.compute_worst_case_cost_controller(
            autom,
            [3];
            initial_set = [1],
            safe_set = [1, 3, 4, 5],            # state 2 is not safe
        )

    @test !(2 in controllable)
    @test 2 in uncontrollable
    @test 1 in controllable
    @test value[1] == 3                     # one step longer: 1 → 4 → 5 → 3
    @test ST.output_control(controller, nothing, 1) == 2   # takes u2, the long way
end

# 1 --u1--> {2, 3}   nondeterministic: the input *may* land in 3
# 3 --u1--> 2        so without an avoid part, going through 3 is fine
# 1 --u2--> 4 --u1--> 2
function leaky_automaton()
    autom = SY.IndexedAutomatonList(4, 2)
    SY.add_transition!(autom, 1, 2, 1)
    SY.add_transition!(autom, 1, 3, 1)
    SY.add_transition!(autom, 3, 2, 1)
    SY.add_transition!(autom, 1, 4, 2)
    SY.add_transition!(autom, 4, 2, 1)
    return autom
end

@testset "an input that may leak into an unsafe state is rejected" begin
    autom = leaky_automaton()
    @test !SY.is_deterministic(autom)

    _, controllable, _, _ =
        OPDS.compute_worst_case_cost_controller(autom, [2]; initial_set = [1])
    @test 3 in controllable                 # nothing forbids passing through 3

    controller, controllable, _, _ = OPDS.compute_worst_case_cost_controller(
        autom,
        [2];
        initial_set = [1],
        safe_set = [1, 2, 4],               # state 3 is not safe
    )

    @test !(3 in controllable)
    @test 1 in controllable
    # `u1` may land in the unsafe state 3, so only `u2` — which cannot — is admissible.
    @test ST.output_control(controller, nothing, 1) == 2
end

@testset "the safe set is honoured on the cost-function path too" begin
    # A non-`nothing` cost routes through `compute_optimal_controller`, a different fixed
    # point (Dijkstra) from the uniform-cost one exercised above.
    autom = detour_automaton()
    unit_cost = (q, u) -> 1.0

    _, _, _, value = OPDS.compute_worst_case_cost_controller(
        autom,
        [3];
        initial_set = [1],
        cost_function = unit_cost,
    )
    @test value[1] == 2

    _, controllable, _, value = OPDS.compute_worst_case_cost_controller(
        autom,
        [3];
        initial_set = [1],
        safe_set = [1, 3, 4, 5],
        cost_function = unit_cost,
    )
    @test !(2 in controllable)
    @test value[1] == 3
end

@testset "a target outside the safe set does not count as reached" begin
    # `safe U target`: the safe part must still hold when the target is entered, so a target
    # state that is not safe is no target at all.
    autom = detour_automaton()
    _, controllable, _, _ = OPDS.compute_worst_case_cost_controller(
        autom,
        [2];
        initial_set = [1],
        safe_set = [1, 3, 4, 5],
    )

    @test isempty(controllable)
end

@testset "the problem carries the safe set through the solver" begin
    autom = detour_automaton()
    problem = PR.OptimalControlProblem(
        autom,
        [1],
        [3],
        nothing,
        nothing,
        PR.Infinity();
        safe_set = [1, 3, 4, 5],
    )
    @test problem.safe_set == [1, 3, 4, 5]

    optimizer = MOI.instantiate(OPDS.OptimizerOptimalControlProblem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    @test optimizer.success
    @test !(2 in optimizer.controllable_set)
    @test optimizer.value_fun_tab[1] == 3
end

@testset "a safe set the solver cannot honour is refused, not ignored" begin
    autom = detour_automaton()
    problem = PR.OptimalControlProblem(
        autom,
        [1],
        [3],
        nothing,
        nothing,
        PR.Infinity();
        safe_set = [1, 3, 4, 5],
    )
    @test_throws ErrorException PR.check_safe_set_supported(problem, "SomeSolver")

    plain = PR.OptimalControlProblem(autom, [1], [3], nothing, nothing, PR.Infinity())
    @test plain.safe_set === nothing
    @test PR.check_safe_set_supported(plain, "SomeSolver") === nothing
end

end # module TestMain
