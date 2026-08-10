module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))
using MathOptInterface
const MOI = MathOptInterface

import MathematicalSystems as MS
import LazySets
using HybridSystems

@testset "HybridSystemAbstraction - simple example, time taken into account" begin
    # ------------------------------
    # Define the HybridSystem
    # ------------------------------

    # Define state and input sets
    X = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    U = LazySets.Hyperrectangle(; low = [-1.5], high = [1.5])

    # Define system dynamics for two modes
    mode1_f(x, u) = [0.5 * x[1] + u[1]]
    mode2_f(x, u) = [0.8 * x[1] + u[1]]

    mode1_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system
    time_sys = MS.ConstrainedLinearContinuousSystem(
        [1.0;;],
        LazySets.Hyperrectangle(; low = [0.0], high = [3.0]),
    )

    guard_1 = LazySets.Hyperrectangle(; low = [0.2, 0.0], high = [1.0, 2.0])

    reset_map = ST.GuardedResetMap(guard_1, _ -> [0.0, 0.0])

    # Automaton and hybrid system
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    modes_systems = [
        ST.VectorContinuousSystem([mode1_system, time_sys]),
        ST.VectorContinuousSystem([mode2_system, time_sys]),
    ]
    reset_maps = [reset_map]
    switchings = [HybridSystems.AutonomousSwitching()]
    concrete_system =
        HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # ------------------------------
    # Define the Control Problem
    # ------------------------------

    initial_state = ([0.0], 0.0, 1) # (state, time, mode)
    Xs_target = [LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))]
    Ts_target = [LazySets.Hyperrectangle(; low = SVector(1.0), high = SVector(2.0))]
    Ns_target = [2]
    target_set = PR.hybrid_reach_spec(Xs_target, Ts_target, Ns_target)
    transition_cost = (aug_state, u) -> 1.0
    concrete_problem = PR.OptimalControlProblem(
        concrete_system,
        initial_state,
        target_set,
        nothing,
        transition_cost,
        PR.Infinity(),
    )

    # ------------------------------
    # Define the Solver parameters
    # ------------------------------

    # Abstraction subsolver parameters (by mode)
    optimizer_list = [
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
    ]

    state_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
    input_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
    state_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))
    input_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.1,
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.1,
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
    ]

    optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
        optimizer_kwargs_dict,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

    # Solve using optimizer
    MOI.optimize!(optimizer)

    # Retrieve results
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # Validate concrete problem
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.transition_cost == transition_cost

    # Validate abstract target set
    abstract_target_set = SY.states_satisfying(abstract_system, target_set)
    for q in abstract_target_set
        (x, t, k) = SY.get_concrete_state(abstract_system, q)
        idx = findfirst(==(k), Ns_target)
        @test !isnothing(idx)
        @test x ∈ Xs_target[idx]
        @test t ≥ LazySets.low(Ts_target[idx], 1) && t ≤ LazySets.high(Ts_target[idx], 1)
    end

    # Validate abstract problem
    @test abstract_problem.initial_set ==
          [SY.get_abstract_state(abstract_system, concrete_problem.initial_set)]
    @test abstract_problem.target_set == abstract_target_set
    @test abstract_problem.state_cost == concrete_problem.state_cost
    @test abstract_problem.time == concrete_problem.time

    for state in 1:SY.get_n_state(abstract_system)
        for input in 1:(abstract_system.input_mapping.total_inputs)
            @test abstract_problem.transition_cost(state, input) == 1.0
        end
    end

    # Validate controllers
    @test !isnothing(abstract_controller)
    abstract_controllable_set = abstract_problem.target_set # Simplified assumption for test
    @test !isnothing(concrete_controller)

    for q in abstract_controllable_set
        (x, t, k) = SY.get_concrete_state(abstract_system, q)
        aug_state = (x, t, k)
        idx = findfirst(==(k), Ns_target)
        in_target =
            !isnothing(idx) &&
            (x ∈ Xs_target[idx]) &&
            (t ≥ LazySets.low(Ts_target[idx], 1)) &&
            (t ≤ LazySets.high(Ts_target[idx], 1))
        if !in_target
            @test concrete_controller.f(aug_state) !== nothing
        end
    end

    # Test reached function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.5, 1.5, 2))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.5, 1.5, 1))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(2.0, 1.5, 2))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.0, 0.5, 2))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.0, 2.5, 2))
    @test AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(-1.0, 1.0, 2))
    @test AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(1.0, 2.0, 2))

    # Test get_next_aug_state
    aug_state = ([0.0], 0.0, 1)
    u_cont = [0.5]
    k = aug_state[3]
    map_sys = ST.simulate_control_map(HybridSystems.mode(concrete_system, k).systems[1].f)
    next_aug_state = AB.HybridSystemAbstraction.get_next_aug_state(
        concrete_system,
        aug_state,
        u_cont,
        true,
        0.1,
        map_sys,
    )

    @test length(next_aug_state) == 3
    u_switch = "SWITCH 1 -> 2"
    next_aug_state_switch = AB.HybridSystemAbstraction.get_next_aug_state(
        concrete_system,
        aug_state,
        u_switch,
        true,
        0.1,
        map_sys,
    )
    @test next_aug_state_switch[3] == 2

    # Test closed-loop trajectory
    reached(aug_x) = AB.HybridSystemAbstraction.reached(concrete_problem, aug_x)

    max_steps = 20
    tsteps = [0.1, 0.1]
    aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
        concrete_system,
        concrete_controller,
        tsteps,
        initial_state,
        max_steps,
        stopping = reached,
    )
    @test !isempty(aug_x_traj)
    @test aug_x_traj[1] == initial_state
    @test length(aug_x_traj) ≤ 21
    @test length(u_traj) ≤ 20
    @test reached(aug_x_traj[end])
    @test length(aug_x_traj) == length(u_traj) + 1
    @test all(x -> length(x) == 3, aug_x_traj)
    @test all(!isnothing, u_traj)
end

@testset "HybridSystemAbstraction - time not taken into account" begin
    # ------------------------------
    # Define the HybridSystem
    # ------------------------------

    # Define state and input sets
    X = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    U = LazySets.Hyperrectangle(; low = [-1.5], high = [1.5])

    # Define system dynamics for two modes
    mode1_f(x, u) = [0.5 * x[1] + u[1]]
    mode2_f(x, u) = [0.8 * x[1] + u[1]]

    mode1_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system (no time evolution)
    time_sys = MS.ConstrainedLinearContinuousSystem(
        [0.0;;],
        LazySets.Hyperrectangle(; low = [0.0], high = [3.0]),
    )

    guard_1 = LazySets.Hyperrectangle(; low = [0.2, 0.0], high = [1.0, 2.0])

    reset_map = ST.GuardedResetMap(guard_1, _ -> [0.0, 0.0])

    # Automaton and hybrid system
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    modes_systems = [
        ST.VectorContinuousSystem([mode1_system, time_sys]),
        ST.VectorContinuousSystem([mode2_system, time_sys]),
    ]
    reset_maps = [reset_map]
    switchings = [HybridSystems.AutonomousSwitching()]
    concrete_system =
        HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # ------------------------------
    # Define the Control Problem
    # ------------------------------

    initial_state = ([0.0], 0.0, 1) # (state, time, mode)
    Xs_target = [LazySets.Hyperrectangle(; low = SVector(-1.0), high = SVector(1.0))]
    Ts_target = [LazySets.Hyperrectangle(; low = SVector(0.0), high = SVector(3.0))]
    Ns_target = [2]
    target_set = PR.hybrid_reach_spec(Xs_target, Ts_target, Ns_target)
    transition_cost = (aug_state, u) -> 1.0
    concrete_problem = PR.OptimalControlProblem(
        concrete_system,
        initial_state,
        target_set,
        nothing,
        transition_cost,
        PR.Infinity(),
    )

    # ------------------------------
    # Define the Solver parameters
    # ------------------------------

    # Abstraction subsolver parameters (by mode)
    optimizer_list = [
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
    ]

    state_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
    input_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
    state_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))
    input_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.1,
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.1,
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
    ]

    optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
        optimizer_kwargs_dict,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

    # Solve using optimizer
    MOI.optimize!(optimizer)

    # Retrieve results
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # Validate concrete problem
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.transition_cost == transition_cost

    # Validate abstract target set
    abstract_target_set = SY.states_satisfying(abstract_system, target_set)
    for q in abstract_target_set
        (x, t, k) = SY.get_concrete_state(abstract_system, q)
        idx = findfirst(==(k), Ns_target)
        @test !isnothing(idx)
        @test x ∈ Xs_target[idx]
        @test t ≥ LazySets.low(Ts_target[idx], 1) && t ≤ LazySets.high(Ts_target[idx], 1)
    end

    # Validate abstract problem
    @test abstract_problem.initial_set ==
          [SY.get_abstract_state(abstract_system, concrete_problem.initial_set)]
    @test abstract_problem.target_set == abstract_target_set
    @test abstract_problem.state_cost == concrete_problem.state_cost
    @test abstract_problem.time == concrete_problem.time

    for state in 1:SY.get_n_state(abstract_system)
        for input in 1:(abstract_system.input_mapping.total_inputs)
            @test abstract_problem.transition_cost(state, input) == 1.0
        end
    end

    # Validate controllers
    @test !isnothing(abstract_controller)
    abstract_controllable_set = abstract_problem.target_set # Simplified assumption for test
    @test !isnothing(concrete_controller)

    for q in abstract_controllable_set
        (x, t, k) = SY.get_concrete_state(abstract_system, q)
        aug_state = (x, t, k)
        idx = findfirst(==(k), Ns_target)
        in_target =
            !isnothing(idx) &&
            (x ∈ Xs_target[idx]) &&
            (t ≥ LazySets.low(Ts_target[idx], 1)) &&
            (t ≤ LazySets.high(Ts_target[idx], 1))
        if !in_target
            @test concrete_controller.h(aug_state) !== nothing
        end
    end

    # Test reached function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.5, 0.0, 2))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.5, 0.0, 1))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(2.0, 0.0, 2))
    @test AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.0, 0.5, 2))
    @test AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(0.0, 2.5, 2))
    @test AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(-1.0, 1.0, 2))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, make_aug_state(1.1, 2.0, 2))

    # Test get_next_aug_state
    aug_state = ([0.0], 0.0, 1)
    u_cont = [0.5]
    k = aug_state[3]
    map_sys = ST.simulate_control_map(HybridSystems.mode(concrete_system, k).systems[1].f)
    next_aug_state = AB.HybridSystemAbstraction.get_next_aug_state(
        concrete_system,
        aug_state,
        u_cont,
        false,
        0.1,
        map_sys,
    )
    @test length(next_aug_state) == 3
    u_switch = "SWITCH 1 -> 2"
    next_aug_state_switch = AB.HybridSystemAbstraction.get_next_aug_state(
        concrete_system,
        aug_state,
        u_switch,
        false,
        0.1,
        map_sys,
    )
    @test next_aug_state_switch[3] == 2

    # Test closed-loop trajectory
    reached(aug_x) = AB.HybridSystemAbstraction.reached(concrete_problem, aug_x)

    max_steps = 20
    tsteps = [0.1, 0.1]
    aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
        concrete_system,
        concrete_controller,
        tsteps,
        initial_state,
        max_steps,
        stopping = reached,
    )

    @test !isempty(aug_x_traj)
    @test aug_x_traj[1] == initial_state
    @test length(aug_x_traj) ≤ 5
    @test length(u_traj) ≤ 4
    @test AB.HybridSystemAbstraction.reached(concrete_problem, aug_x_traj[end])
    @test length(aug_x_traj) == length(u_traj) + 1
    @test all(x -> length(x) == 3, aug_x_traj)
    @test all(!isnothing, u_traj)
end

@testset "HybridSystemAbstraction - safety problem" begin
    # ------------------------------
    # Define the HybridSystem
    # ------------------------------

    # Define a simple 1D system
    X = LazySets.Hyperrectangle(; low = [-10.0], high = [10.0])
    U = LazySets.Hyperrectangle(; low = [-5.5], high = [8.5])

    # Define dynamics
    mode1_f(x, u) = [-0.1 * x[1] + u[1]]
    mode2_f(x, u) = [0.4 * x[1] + u[1]]

    mode1_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system
    time_sys = MS.ConstrainedLinearContinuousSystem(
        [1.0;;],
        LazySets.Hyperrectangle(; low = [0.0], high = [5.0]),
    )

    guard_1 = LazySets.Hyperrectangle(; low = [0.0, 0.0], high = [10.0, 5.0])
    reset_map = ST.GuardedResetMap(guard_1)      # identity on the augmented [x; t]

    # Automaton and hybrid system
    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    modes_systems = [
        ST.VectorContinuousSystem([mode1_system, time_sys]),
        ST.VectorContinuousSystem([mode2_system, time_sys]),
    ]
    reset_maps = [reset_map]
    switchings = [HybridSystems.AutonomousSwitching()]
    concrete_system =
        HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # ------------------------------
    # Define the Control Problem
    # ------------------------------

    initial_state = ([0.0], 0.0, 1) # (state, time, mode)
    Xs_safe = [
        LazySets.Hyperrectangle(; low = [-3.0], high = [3.0]),
        LazySets.Hyperrectangle(; low = [-1.0], high = [10.0]),
    ]
    Ts_safe = [
        LazySets.Hyperrectangle(; low = [0.0], high = [2.0]),
        LazySets.Hyperrectangle(; low = [1.0], high = [5.0]),
    ]
    Ns_safe = [1, 2]
    safe_set = PR.hybrid_reach_spec(Xs_safe, Ts_safe, Ns_safe; incl_mode = UT.OUTER)
    concrete_problem = PR.SafetyProblem(concrete_system, initial_state, safe_set, 10.0)

    # ------------------------------
    # Define the Solver parameters
    # ------------------------------

    # Abstraction parameters
    optimizer_list = [
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
    ]

    state_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
    input_grid_1 = MP.GridFree(SVector(0.0), SVector(0.1))
    state_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))
    input_grid_2 = MP.GridFree(SVector(0.0), SVector(0.1))

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.1,
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.1),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.1,
            "approx_mode" => AB.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.4),
        ),
    ]

    optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
        optimizer_kwargs_dict,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

    # Solve using optimizer
    MOI.optimize!(optimizer)

    # Retrieve results
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # Validate concrete problem
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.time == 10.0
    @test isa(concrete_problem, PR.SafetyProblem)

    # Validate abstract safe set
    abstract_safe_set = SY.states_satisfying(abstract_system, safe_set)
    for q in abstract_safe_set
        (x, t, k) = SY.get_concrete_state(abstract_system, q)
        idx = findfirst(==(k), Ns_safe)
        @test !isnothing(idx)
        @test x ∈ Xs_safe[idx]
        @test t ≥ LazySets.low(Ts_safe[idx], 1) && t ≤ LazySets.high(Ts_safe[idx], 1)
    end

    # Validate abstract problem
    @test abstract_problem.initial_set ==
          [SY.get_abstract_state(abstract_system, concrete_problem.initial_set)]
    @test abstract_problem.safe_set == abstract_safe_set
    @test abstract_problem.time == concrete_problem.time
    @test isa(abstract_problem, PR.SafetyProblem)

    # Validate controllers
    @test !isnothing(abstract_controller)
    @test !isnothing(concrete_controller)

    # Test safety checking function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.HybridSystemAbstraction.safe(concrete_problem, make_aug_state(0.0, 1.0, 1))
    @test AB.HybridSystemAbstraction.safe(concrete_problem, make_aug_state(2.0, 3.0, 2))
    @test !AB.HybridSystemAbstraction.safe(concrete_problem, make_aug_state(4.0, 1.0, 1))
    @test AB.HybridSystemAbstraction.safe(concrete_problem, make_aug_state(-0.9, 2.0, 2))

    # Test closed-loop trajectory
    unsafe(aug_x) = !AB.HybridSystemAbstraction.safe(concrete_problem, aug_x)

    max_steps = 100000
    tsteps = [0.1, 0.1]
    aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
        concrete_system,
        concrete_controller,
        tsteps,
        initial_state,
        max_steps,
        stopping = unsafe,
    )

    @test !isempty(aug_x_traj)
    @test aug_x_traj[1] == initial_state
    @test length(aug_x_traj) == length(u_traj) + 1
    @test all(x -> length(x) == 3, aug_x_traj)
    for state in aug_x_traj[1:(end - 1)]
        @test AB.HybridSystemAbstraction.safe(concrete_problem, state)
    end
    println(
        "Safety test completed: trajectory remained safe for $(length(aug_x_traj)) steps",
    )
end

@testset "HybridSystemAbstraction - time-free hybrid (plain modes)" begin
    # Modes are plain physical systems: NO time subsystem, so the abstraction has no
    # time axis and the augmented state is (x, mode).
    X = LazySets.Hyperrectangle(; low = [-1.0], high = [1.0])
    U = LazySets.Hyperrectangle(; low = [-1.5], high = [1.5])
    mode1_f(x, u) = [0.5 * x[1] + u[1]]
    mode2_f(x, u) = [0.8 * x[1] + u[1]]
    mode1_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system = MS.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    guard_1 = LazySets.Hyperrectangle(; low = [0.2], high = [1.0])            # guard over x only
    reset_map = ST.GuardedResetMap(guard_1, _ -> [0.0])

    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)
    modes_systems = [mode1_system, mode2_system]   # plain, not VectorContinuousSystem
    switchings = [HybridSystems.AutonomousSwitching()]
    concrete_system =
        HybridSystems.HybridSystem(automaton, modes_systems, [reset_map], switchings)

    # (x, mode) augmented state — no time. Target: anywhere in mode 2.
    initial_state = ([0.0], 1)
    target_spec = PR.HybridSpec(
        Dict(2 => PR.StateSpec(LazySets.Hyperrectangle(; low = [-1.0], high = [1.0]))),
    )
    transition_cost = (aug_state, u) -> 1.0
    concrete_problem = PR.OptimalControlProblem(
        concrete_system,
        initial_state,
        target_spec,
        nothing,
        transition_cost,
        PR.Infinity(),
    )

    optimizer_list = [
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(AB.UniformGridAbstraction.Optimizer),
    ]
    grid = MP.GridFree(SVector(0.0), SVector(0.1))
    make_kwargs(jac) = Dict(
        "state_grid" => grid,
        "input_grid" => grid,
        "time_step" => 0.1,                    # integration step for the dynamics
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => u -> SMatrix{1, 1}(jac),
    )
    optimizer_kwargs_dict = [make_kwargs(0.5), make_kwargs(0.8)]

    optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
        optimizer_kwargs_dict,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)

    abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # The abstraction is genuinely time-free: no mode is clock-lifted.
    @test all(m -> m isa SY.SymbolicModelList, abstract_system.mode_models)

    # Augmented states are (x, mode) 2-tuples, and no time factor inflates the count:
    # the flattened state count equals the (local_id, mode) pairs of the plain modes.
    q1 = first(SY.enum_states(abstract_system))
    @test length(SY.get_concrete_state(abstract_system, q1)) == 2

    # Concrete ↔ abstract round-trip on the (x, mode) state.
    for q in SY.enum_states(abstract_system)
        (x, k) = SY.get_concrete_state(abstract_system, q)
        @test SY.get_abstract_state(abstract_system, (x, k)) == q
    end

    @test !isnothing(concrete_controller)

    # Out-of-domain query (x far outside the grid) is handled gracefully, not a crash.
    @test !ST.is_defined(concrete_controller, nothing, ([100.0], 1))
    @test ST.output_control(concrete_controller, nothing, ([100.0], 1)) === nothing

    # reached uses the arity-agnostic spec membership.
    @test AB.HybridSystemAbstraction.reached(concrete_problem, ([0.0], 2))
    @test !AB.HybridSystemAbstraction.reached(concrete_problem, ([0.0], 1))

    # Time-free closed-loop: the augmented state stays a (x, mode) 2-tuple throughout.
    reached_fn(aug) = AB.HybridSystemAbstraction.reached(concrete_problem, aug)
    aug_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
        concrete_system,
        concrete_controller,
        [0.1, 0.1],
        initial_state,
        50;
        stopping = reached_fn,
    )
    @test aug_traj[1] == initial_state
    @test all(s -> length(s) == 2, aug_traj)
    @test length(aug_traj) == length(u_traj) + 1
    @test reached_fn(aug_traj[end])
end

end # module
