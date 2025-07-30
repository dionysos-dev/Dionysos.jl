module TestMain
using Test     #src

using StaticArrays, Plots, HybridSystems, MathematicalSystems

# At this point, we import the useful Dionysos sub-modules.
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

@testset "Timedhybrid_abstraction - simple exemple, time taken into account" begin
    # Define state and input sets
    X = UT.HyperRectangle([-1.0], [1.0])
    U = UT.HyperRectangle([-1.5], [1.5])

    # Define system dynamics for two modes
    mode1_f(x, u) = [0.5 * x[1] + u[1]]
    mode2_f(x, u) = [0.8 * x[1] + u[1]]

    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system
    time_sys = MathematicalSystems.ConstrainedLinearContinuousSystem([1.0;;], UT.HyperRectangle([0.0], [3.0]))

    # Reset map
    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end
    MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector) = reset.target
    MathematicalSystems.stateset(reset::FixedPointResetMap) = reset.domain

    guard_1 = UT.HyperRectangle([0.2, 0.0], [1.0, 2.0])
    reset_map = FixedPointResetMap(guard_1, [0.0, 0.0])

    # Automaton and hybrid system
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    modes_systems = [
        SY.VectorContinuousSystem([mode1_system, time_sys]),
        SY.VectorContinuousSystem([mode2_system, time_sys]),
    ]
    reset_maps = [reset_map]
    switchings = [HybridSystems.AutonomousSwitching()]
    hs = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Abstraction parameters
    growth_bounds = SVector(SMatrix{1, 1}(0.5), SMatrix{1, 1}(0.8))
    param_discretization = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.1)]
    hybrid_symmodel = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(hs, growth_bounds, param_discretization)

    # Problem specification
    initial_state = ([0.0], 0.0, 1)
    Xs_target = [UT.HyperRectangle([-1.0], [1.0])]
    Ts_target = [UT.HyperRectangle([1.0], [2.0])]
    Ns_target = [2]
    cost_fun = (aug_state, u) -> 1.0
    concret_specs = AB.TemporalHybridSymbolicModelAbstraction.ProblemSpecs(initial_state, Xs_target, Ts_target, Ns_target, cost_fun)

    # Concrete problem
    concrete_problem = AB.TemporalHybridSymbolicModelAbstraction.build_concrete_problem(concret_specs)
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.transition_cost == cost_fun
    @test concrete_problem.time == Dionysos.Problem.Infinity()

    # Abstract target set
    abstract_target_set = SY.TimedHybridAutomata.get_states_from_set(hybrid_symmodel, Xs_target, Ts_target, Ns_target)
    for q in abstract_target_set
        (x, t, k) = SY.TimedHybridAutomata.get_concrete_state(hybrid_symmodel, q)
        idx = findfirst(==(k), Ns_target)
        @test !isnothing(idx)
        @test x ∈ Xs_target[idx]
        @test t ≥ Ts_target[idx].lb[1] && t ≤ Ts_target[idx].ub[1]
    end

    # Abstract problem
    abstract_problem = AB.TemporalHybridSymbolicModelAbstraction.build_abstract_problem(concrete_problem, hybrid_symmodel)
    @test abstract_problem.initial_set == [SY.TimedHybridAutomata.get_abstract_state(hybrid_symmodel, concrete_problem.initial_set)]
    @test abstract_problem.target_set == abstract_target_set
    @test abstract_problem.state_cost == concrete_problem.state_cost
    @test abstract_problem.time == concrete_problem.time

    for state in 1:SY.TimedHybridAutomata.get_n_state(hybrid_symmodel)
        for input in 1:hybrid_symmodel.global_input_map.total_inputs
            @test abstract_problem.transition_cost(state, input) == 1.0
        end
    end

    # Solve abstract and concrete problems
    abstract_controller, controllable_set_symbols = AB.TemporalHybridSymbolicModelAbstraction.solve_abstract_problem(abstract_problem)
    @test !isnothing(abstract_controller)
    @test !isempty(controllable_set_symbols)

    concrete_controller = AB.TemporalHybridSymbolicModelAbstraction.solve_concrete_problem(hybrid_symmodel, abstract_controller)
    for q in controllable_set_symbols
        (x, t, k) = SY.TimedHybridAutomata.get_concrete_state(hybrid_symmodel, q)
        aug_state = (x, t, k)
        idx = findfirst(==(k), Ns_target)
        in_target = !isnothing(idx) && (x ∈ Xs_target[idx]) && (t ≥ Ts_target[idx].lb[1]) && (t ≤ Ts_target[idx].ub[1])
        if !in_target
            @test concrete_controller.f(aug_state) !== nothing
        end
    end

    # Test reached function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.5, 1.5, 2))
    @test !AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.5, 1.5, 1))
    @test !AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(2.0, 1.5, 2))
    @test !AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.0, 0.5, 2))
    @test !AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.0, 2.5, 2))
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(-1.0, 1.0, 2))
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(1.0, 2.0, 2))

    # Test get_next_aug_state
    aug_state = ([0.0], 0.0, 1)
    u_cont = [0.5]
    k = aug_state[3]
    tm = hybrid_symmodel.time_symbolic_models[k]
    map_sys = ST.simulate_control_map(HybridSystems.mode(hs, k).systems[1].f)
    next_aug_state = AB.TemporalHybridSymbolicModelAbstraction.get_next_aug_state(hs, aug_state, u_cont, tm, map_sys)
    @test length(next_aug_state) == 3
    u_switch = "SWITCH 1 -> 2"
    next_aug_state_switch = AB.TemporalHybridSymbolicModelAbstraction.get_next_aug_state(hs, aug_state, u_switch, tm, map_sys)
    @test next_aug_state_switch[3] == 2

    # Test closed-loop trajectory
    traj, ctrls = AB.TemporalHybridSymbolicModelAbstraction.get_closed_loop_trajectory(
        hybrid_symmodel, hs, concret_specs, concrete_controller, initial_state, 20,
        stopping = AB.TemporalHybridSymbolicModelAbstraction.reached,
    )
    @test !isempty(traj)
    @test traj[1] == initial_state
    @test length(traj) ≤ 21
    @test length(ctrls) ≤ 20
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, traj[end])
    @test length(traj) == length(ctrls) + 1
    @test all(x -> length(x) == 3, traj)
    @test all(!isnothing, ctrls)

    # Test solve shortcut
    controller = AB.TemporalHybridSymbolicModelAbstraction.solve(
        hs, growth_bounds, param_discretization, concret_specs,
    )
    for state in traj[1:end-1]
        @test controller.f(state) == concrete_controller.f(state)
    end
end

@testset "Timedhybrid_abstraction - time not taken into account" begin
    # Define state and input sets
    X = UT.HyperRectangle([-1.0], [1.0])
    U = UT.HyperRectangle([-1.5], [1.5])

    # Define system dynamics for two modes
    mode1_f(x, u) = [0.5 * x[1] + u[1]]
    mode2_f(x, u) = [0.8 * x[1] + u[1]]

    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system
    time_sys = MathematicalSystems.ConstrainedLinearContinuousSystem([0.0;;], UT.HyperRectangle([0.0], [3.0]))

    # Reset map
    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end
    MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector) = reset.target
    MathematicalSystems.stateset(reset::FixedPointResetMap) = reset.domain

    guard_1 = UT.HyperRectangle([0.2, 0.0], [1.0, 2.0])
    reset_map = FixedPointResetMap(guard_1, [0.0, 0.0])

    # Automaton and hybrid system
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    modes_systems = [
        SY.VectorContinuousSystem([mode1_system, time_sys]),
        SY.VectorContinuousSystem([mode2_system, time_sys]),
    ]
    reset_maps = [reset_map]
    switchings = [HybridSystems.AutonomousSwitching()]
    hs = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Abstraction parameters
    growth_bounds = SVector(SMatrix{1, 1}(0.5), SMatrix{1, 1}(0.8))
    param_discretization = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.1)]
    hybrid_symmodel = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(hs, growth_bounds, param_discretization)

    # Problem specification
    initial_state = ([0.0], 0.0, 1)
    Xs_target = [UT.HyperRectangle([-1.0], [1.0])]
    Ts_target = [UT.HyperRectangle([0.0], [3.0])] # Time is not taken into account
    Ns_target = [2]
    cost_fun = (aug_state, u) -> 1.0
    concret_specs = AB.TemporalHybridSymbolicModelAbstraction.ProblemSpecs(initial_state, Xs_target, Ts_target, Ns_target, cost_fun)

    # Concrete problem
    concrete_problem = AB.TemporalHybridSymbolicModelAbstraction.build_concrete_problem(concret_specs)
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.transition_cost == cost_fun
    @test concrete_problem.time == Dionysos.Problem.Infinity()

    # Abstract target set
    abstract_target_set = SY.TimedHybridAutomata.get_states_from_set(hybrid_symmodel, Xs_target, Ts_target, Ns_target)
    for q in abstract_target_set
        (x, t, k) = SY.TimedHybridAutomata.get_concrete_state(hybrid_symmodel, q)
        idx = findfirst(==(k), Ns_target)
        @test !isnothing(idx)
        @test x ∈ Xs_target[idx]
        @test t ≥ Ts_target[idx].lb[1] && t ≤ Ts_target[idx].ub[1]
    end

    # Abstract problem
    abstract_problem = AB.TemporalHybridSymbolicModelAbstraction.build_abstract_problem(concrete_problem, hybrid_symmodel)
    @test abstract_problem.initial_set == [SY.TimedHybridAutomata.get_abstract_state(hybrid_symmodel, concrete_problem.initial_set)]
    @test abstract_problem.target_set == abstract_target_set
    @test abstract_problem.state_cost == concrete_problem.state_cost
    @test abstract_problem.time == concrete_problem.time

    for state in 1:SY.TimedHybridAutomata.get_n_state(hybrid_symmodel)
        for input in 1:hybrid_symmodel.global_input_map.total_inputs
            @test abstract_problem.transition_cost(state, input) == 1.0
        end
    end

    # Solve abstract and concrete problems
    abstract_controller, controllable_set_symbols = AB.TemporalHybridSymbolicModelAbstraction.solve_abstract_problem(abstract_problem)
    @test !isnothing(abstract_controller)
    @test !isempty(controllable_set_symbols)

    concrete_controller = AB.TemporalHybridSymbolicModelAbstraction.solve_concrete_problem(hybrid_symmodel, abstract_controller)
    for q in controllable_set_symbols
        (x, t, k) = SY.TimedHybridAutomata.get_concrete_state(hybrid_symmodel, q)
        aug_state = (x, t, k)
        idx = findfirst(==(k), Ns_target)
        in_target = !isnothing(idx) && (x ∈ Xs_target[idx]) && (t ≥ Ts_target[idx].lb[1]) && (t ≤ Ts_target[idx].ub[1])
        if !in_target
            @test concrete_controller.f(aug_state) !== nothing
        end
    end

    # Test reached function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.5, 0.0, 2))
    @test !AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.5, 0.0, 1))
    @test !AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(2.0, 0.0, 2))
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.0, 0.5, 2))
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(0.0, 2.5, 2))
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(-1.0, 1.0, 2))
    @test !AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, make_aug_state(1.1, 2.0, 2))

    # Test get_next_aug_state
    aug_state = ([0.0], 0.0, 1)
    u_cont = [0.5]
    k = aug_state[3]
    tm = hybrid_symmodel.time_symbolic_models[k]
    map_sys = ST.simulate_control_map(HybridSystems.mode(hs, k).systems[1].f)
    next_aug_state = AB.TemporalHybridSymbolicModelAbstraction.get_next_aug_state(hs, aug_state, u_cont, tm, map_sys)
    @test length(next_aug_state) == 3
    u_switch = "SWITCH 1 -> 2"
    next_aug_state_switch = AB.TemporalHybridSymbolicModelAbstraction.get_next_aug_state(hs, aug_state, u_switch, tm, map_sys)
    @test next_aug_state_switch[3] == 2

    # Test closed-loop trajectory
    traj, ctrls = AB.TemporalHybridSymbolicModelAbstraction.get_closed_loop_trajectory(
        hybrid_symmodel, hs, concret_specs, concrete_controller, initial_state, 20,
        stopping = AB.TemporalHybridSymbolicModelAbstraction.reached,
    )
    
    println("Trajectory: ", traj)
    println("Controllers: ", ctrls)

    @test !isempty(traj)
    @test traj[1] == initial_state
    @test length(traj) ≤ 4
    @test length(ctrls) ≤ 3
    @test AB.TemporalHybridSymbolicModelAbstraction.reached(concret_specs, traj[end])
    @test length(traj) == length(ctrls) + 1
    @test all(x -> length(x) == 3, traj)
    @test all(!isnothing, ctrls)

    # Test solve shortcut
    controller = AB.TemporalHybridSymbolicModelAbstraction.solve(
        hs, growth_bounds, param_discretization, concret_specs,
    )
    for state in traj[1:end-1]
        @test controller.f(state) == concrete_controller.f(state)
    end
end
end