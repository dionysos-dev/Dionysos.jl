module TestMain
using Test
using StaticArrays, HybridSystems, MathematicalSystems
using Dionysos
using MathOptInterface

const DI = Dionysos
const MOI = MathOptInterface
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

@testset "Timedhybrid_abstraction - simple example, time taken into account" begin
    # Define state and input sets
    X = UT.HyperRectangle([-1.0], [1.0])
    U = UT.HyperRectangle([-1.5], [1.5])

    # Define system dynamics for two modes
    mode1_f(x, u) = [0.5 * x[1] + u[1]]
    mode2_f(x, u) = [0.8 * x[1] + u[1]]

    mode1_system =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system
    time_sys = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([0.0], [3.0]),
    )

    # Reset map
    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end
    MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector) =
        reset.target
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
    optimizer_factory_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]

    state_grid_1 = DO.GridFree(SVector(0.0), SVector(0.1))
    input_grid_1 = DO.GridFree(SVector(0.0), SVector(0.1))
    state_grid_2 = DO.GridFree(SVector(0.0), SVector(0.1))
    input_grid_2 = DO.GridFree(SVector(0.0), SVector(0.1))

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.1,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.1,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
    ]

    param_discretization = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.1)]

    # Problem specification
    initial_state = ([0.0], 0.0, 1)
    Xs_target = [UT.HyperRectangle([-1.0], [1.0])]
    Ts_target = [UT.HyperRectangle([1.0], [2.0])]
    Ns_target = [2]
    cost_fun = (aug_state, u) -> 1.0
    concret_specs = AB.TimedHybridAbstraction.TimedHybridOptimalControlProblem(
        initial_state,
        Xs_target,
        Ts_target,
        Ns_target,
        cost_fun,
        Dionysos.Problem.Infinity(),
    )

    # Create and configure optimizer
    optimizer = AB.TimedHybridAbstraction.Optimizer()
    AB.TimedHybridAbstraction.set_optimizer!(
        optimizer,
        concret_specs,
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
    )

    # Solve using optimizer
    MOI.optimize!(optimizer)

    # Retrieve results
    concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"))
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    hybrid_symmodel = MOI.get(optimizer, MOI.RawOptimizerAttribute("symbolic_model"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # Validate concrete problem
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.transition_cost == cost_fun

    # Validate abstract target set
    abstract_target_set = SY.SymbolicTimedHybridSystems.get_states_from_set(
        hybrid_symmodel,
        Xs_target,
        Ts_target,
        Ns_target,
    )
    for q in abstract_target_set
        (x, t, k) = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_symmodel, q)
        idx = findfirst(==(k), Ns_target)
        @test !isnothing(idx)
        @test x ∈ Xs_target[idx]
        @test t ≥ Ts_target[idx].lb[1] && t ≤ Ts_target[idx].ub[1]
    end

    # Validate abstract problem
    @test abstract_problem.initial_set == [
        SY.SymbolicTimedHybridSystems.get_abstract_state(
            hybrid_symmodel,
            concrete_problem.initial_set,
        ),
    ]
    @test abstract_problem.target_set == abstract_target_set
    @test abstract_problem.state_cost == concrete_problem.state_cost
    @test abstract_problem.time == concrete_problem.time

    for state in 1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_symmodel)
        for input in 1:hybrid_symmodel.input_mapping.total_inputs
            @test abstract_problem.transition_cost(state, input) == 1.0
        end
    end

    # Validate controllers
    @test !isnothing(abstract_controller)
    controllable_set_symbols = abstract_problem.target_set # Simplified assumption for test
    @test !isnothing(concrete_controller)

    for q in controllable_set_symbols
        (x, t, k) = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_symmodel, q)
        aug_state = (x, t, k)
        idx = findfirst(==(k), Ns_target)
        in_target =
            !isnothing(idx) &&
            (x ∈ Xs_target[idx]) &&
            (t ≥ Ts_target[idx].lb[1]) &&
            (t ≤ Ts_target[idx].ub[1])
        if !in_target
            @test concrete_controller.f(aug_state) !== nothing
        end
    end

    # Test reached function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.5, 1.5, 2))
    @test !AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.5, 1.5, 1))
    @test !AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(2.0, 1.5, 2))
    @test !AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.0, 0.5, 2))
    @test !AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.0, 2.5, 2))
    @test AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(-1.0, 1.0, 2))
    @test AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(1.0, 2.0, 2))

    # Test get_next_aug_state
    aug_state = ([0.0], 0.0, 1)
    u_cont = [0.5]
    k = aug_state[3]
    tm = hybrid_symmodel.time_abstractions[k]
    map_sys = ST.simulate_control_map(HybridSystems.mode(hs, k).systems[1].f)
    next_aug_state = AB.TimedHybridAbstraction.get_next_aug_state(
        hs,
        aug_state,
        u_cont,
        true,
        0.1,
        map_sys,
    )
    @test length(next_aug_state) == 3
    u_switch = "SWITCH 1 -> 2"
    next_aug_state_switch = AB.TimedHybridAbstraction.get_next_aug_state(
        hs,
        aug_state,
        u_switch,
        true,
        0.1,
        map_sys,
    )
    @test next_aug_state_switch[3] == 2

    # Test closed-loop trajectory
    traj, ctrls = AB.TimedHybridAbstraction.get_closed_loop_trajectory(
        param_discretization,
        hs,
        concret_specs,
        concrete_controller,
        initial_state,
        20,
        stopping = AB.TimedHybridAbstraction.reached,
    )
    @test !isempty(traj)
    @test traj[1] == initial_state
    @test length(traj) ≤ 21
    @test length(ctrls) ≤ 20
    @test AB.TimedHybridAbstraction.reached(concret_specs, traj[end])
    @test length(traj) == length(ctrls) + 1
    @test all(x -> length(x) == 3, traj)
    @test all(!isnothing, ctrls)

    # Test solve shortcut
    controller = AB.TimedHybridAbstraction.solve_timed_hybrid_problem(
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
        concret_specs,
    )
    for state in traj[1:(end - 1)]
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

    mode1_system =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system (no time evolution)
    time_sys = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [0.0;;],
        UT.HyperRectangle([0.0], [3.0]),
    )

    # Reset map
    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end
    MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector) =
        reset.target
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
    optimizer_factory_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]

    state_grid_1 = DO.GridFree(SVector(0.0), SVector(0.1))
    input_grid_1 = DO.GridFree(SVector(0.0), SVector(0.1))
    state_grid_2 = DO.GridFree(SVector(0.0), SVector(0.1))
    input_grid_2 = DO.GridFree(SVector(0.0), SVector(0.1))

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.1,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.1,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
    ]

    param_discretization = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.1)]

    # Problem specification
    initial_state = ([0.0], 0.0, 1)
    Xs_target = [UT.HyperRectangle([-1.0], [1.0])]
    Ts_target = [UT.HyperRectangle([0.0], [3.0])]
    Ns_target = [2]
    cost_fun = (aug_state, u) -> 1.0
    concret_specs = AB.TimedHybridAbstraction.TimedHybridOptimalControlProblem(
        initial_state,
        Xs_target,
        Ts_target,
        Ns_target,
        cost_fun,
        Dionysos.Problem.Infinity(),
    )

    # Create and configure optimizer
    optimizer = AB.TimedHybridAbstraction.Optimizer()
    AB.TimedHybridAbstraction.set_optimizer!(
        optimizer,
        concret_specs,
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
    )

    # Solve using optimizer
    MOI.optimize!(optimizer)

    # Retrieve results
    concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"))
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    hybrid_symmodel = MOI.get(optimizer, MOI.RawOptimizerAttribute("symbolic_model"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # Validate concrete problem
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.transition_cost == cost_fun

    # Validate abstract target set
    abstract_target_set = SY.SymbolicTimedHybridSystems.get_states_from_set(
        hybrid_symmodel,
        Xs_target,
        Ts_target,
        Ns_target,
    )
    for q in abstract_target_set
        (x, t, k) = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_symmodel, q)
        idx = findfirst(==(k), Ns_target)
        @test !isnothing(idx)
        @test x ∈ Xs_target[idx]
        @test t ≥ Ts_target[idx].lb[1] && t ≤ Ts_target[idx].ub[1]
    end

    # Validate abstract problem
    @test abstract_problem.initial_set == [
        SY.SymbolicTimedHybridSystems.get_abstract_state(
            hybrid_symmodel,
            concrete_problem.initial_set,
        ),
    ]
    @test abstract_problem.target_set == abstract_target_set
    @test abstract_problem.state_cost == concrete_problem.state_cost
    @test abstract_problem.time == concrete_problem.time

    for state in 1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_symmodel)
        for input in 1:hybrid_symmodel.input_mapping.total_inputs
            @test abstract_problem.transition_cost(state, input) == 1.0
        end
    end

    # Validate controllers
    @test !isnothing(abstract_controller)
    controllable_set_symbols = abstract_problem.target_set # Simplified assumption for test
    @test !isnothing(concrete_controller)

    for q in controllable_set_symbols
        (x, t, k) = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_symmodel, q)
        aug_state = (x, t, k)
        idx = findfirst(==(k), Ns_target)
        in_target =
            !isnothing(idx) &&
            (x ∈ Xs_target[idx]) &&
            (t ≥ Ts_target[idx].lb[1]) &&
            (t ≤ Ts_target[idx].ub[1])
        if !in_target
            @test concrete_controller.f(aug_state) !== nothing
        end
    end

    # Test reached function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.5, 0.0, 2))
    @test !AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.5, 0.0, 1))
    @test !AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(2.0, 0.0, 2))
    @test AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.0, 0.5, 2))
    @test AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(0.0, 2.5, 2))
    @test AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(-1.0, 1.0, 2))
    @test !AB.TimedHybridAbstraction.reached(concret_specs, make_aug_state(1.1, 2.0, 2))

    # Test get_next_aug_state
    aug_state = ([0.0], 0.0, 1)
    u_cont = [0.5]
    k = aug_state[3]
    tm = hybrid_symmodel.time_abstractions[k]
    map_sys = ST.simulate_control_map(HybridSystems.mode(hs, k).systems[1].f)
    next_aug_state = AB.TimedHybridAbstraction.get_next_aug_state(
        hs,
        aug_state,
        u_cont,
        false,
        0.1,
        map_sys,
    )
    @test length(next_aug_state) == 3
    u_switch = "SWITCH 1 -> 2"
    next_aug_state_switch = AB.TimedHybridAbstraction.get_next_aug_state(
        hs,
        aug_state,
        u_switch,
        false,
        0.1,
        map_sys,
    )
    @test next_aug_state_switch[3] == 2

    # Test closed-loop trajectory
    traj, ctrls = AB.TimedHybridAbstraction.get_closed_loop_trajectory(
        param_discretization,
        hs,
        concret_specs,
        concrete_controller,
        initial_state,
        20,
        stopping = AB.TimedHybridAbstraction.reached,
    )
    println("Trajectory: ", traj)
    println("Controllers: ", ctrls)
    @test !isempty(traj)
    @test traj[1] == initial_state
    @test length(traj) ≤ 4
    @test length(ctrls) ≤ 3
    @test AB.TimedHybridAbstraction.reached(concret_specs, traj[end])
    @test length(traj) == length(ctrls) + 1
    @test all(x -> length(x) == 3, traj)
    @test all(!isnothing, ctrls)

    # Test solve shortcut
    controller = AB.TimedHybridAbstraction.solve_timed_hybrid_problem(
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
        concret_specs,
    )
    for state in traj[1:(end - 1)]
        @test controller.f(state) == concrete_controller.f(state)
    end
end

@testset "Timedhybrid_abstraction - safety problem" begin
    # Define a simple 1D system
    X = UT.HyperRectangle([-10.0], [10.0])
    U = UT.HyperRectangle([-5.5], [8.5])

    # Define dynamics
    mode1_f(x, u) = [-0.1 * x[1] + u[1]]
    mode2_f(x, u) = [0.4 * x[1] + u[1]]

    mode1_system =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode1_f, 1, 1, X, U)
    mode2_system =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(mode2_f, 1, 1, X, U)

    # Time system
    time_sys = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([0.0], [5.0]),
    )

    # Reset map
    struct SafetyResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end
    MathematicalSystems.apply(reset::SafetyResetMap, state::AbstractVector) =
        [state[1], state[2]]
    MathematicalSystems.stateset(reset::SafetyResetMap) = reset.domain

    guard_1 = UT.HyperRectangle([0.0, 0.0], [10.0, 5.0])
    reset_map = SafetyResetMap(guard_1, [0.0, 4.0])

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
    optimizer_factory_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]

    state_grid_1 = DO.GridFree(SVector(0.0), SVector(0.1))
    input_grid_1 = DO.GridFree(SVector(0.0), SVector(0.1))
    state_grid_2 = DO.GridFree(SVector(0.0), SVector(0.1))
    input_grid_2 = DO.GridFree(SVector(0.0), SVector(0.1))

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.1,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.1),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.1,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.4),
        ),
    ]

    param_discretization = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.1)]

    # Safety problem specification
    initial_state = ([0.0], 0.0, 1)
    Xs_safe = [UT.HyperRectangle([-3.0], [3.0]), UT.HyperRectangle([-1.0], [10.0])]
    Ts_safe = [UT.HyperRectangle([0.0], [2.0]), UT.HyperRectangle([1.0], [5.0])]
    Ns_safe = [1, 2]
    safety_specs = AB.TimedHybridAbstraction.TimedHybridSafetyProblem(
        initial_state,
        Xs_safe,
        Ts_safe,
        Ns_safe,
        10.0,
    )

    # Verify safety problem
    @test safety_specs.problem_type == :safety

    # Create and configure optimizer
    optimizer = AB.TimedHybridAbstraction.Optimizer()
    AB.TimedHybridAbstraction.set_optimizer!(
        optimizer,
        safety_specs,
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
    )

    # Solve using optimizer
    MOI.optimize!(optimizer)

    # Retrieve results
    concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"))
    abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
    hybrid_symmodel = MOI.get(optimizer, MOI.RawOptimizerAttribute("symbolic_model"))
    abstract_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
    concrete_controller =
        MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

    # Validate concrete problem
    @test concrete_problem.initial_set == initial_state
    @test concrete_problem.time == 10.0
    @test isa(concrete_problem, Dionysos.Problem.SafetyProblem)

    # Validate abstract safe set
    abstract_safe_set = SY.SymbolicTimedHybridSystems.get_states_from_set(
        hybrid_symmodel,
        Xs_safe,
        Ts_safe,
        Ns_safe,
        domain = Dionysos.Domain.OUTER,
    )
    for q in abstract_safe_set
        (x, t, k) = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_symmodel, q)
        idx = findfirst(==(k), Ns_safe)
        @test !isnothing(idx)
        @test x ∈ Xs_safe[idx]
        @test t ≥ Ts_safe[idx].lb[1] && t ≤ Ts_safe[idx].ub[1]
    end

    # Validate abstract problem
    @test abstract_problem.initial_set == [
        SY.SymbolicTimedHybridSystems.get_abstract_state(
            hybrid_symmodel,
            concrete_problem.initial_set,
        ),
    ]
    @test abstract_problem.safe_set == abstract_safe_set
    @test abstract_problem.time == concrete_problem.time
    @test isa(abstract_problem, Dionysos.Problem.SafetyProblem)

    # Validate controllers
    @test !isnothing(abstract_controller)
    @test !isnothing(concrete_controller)

    # Test closed-loop trajectory
    traj, ctrls = AB.TimedHybridAbstraction.get_closed_loop_trajectory(
        param_discretization,
        hs,
        safety_specs,
        concrete_controller,
        initial_state,
        100000,
        stopping = (specs, state) -> !AB.TimedHybridAbstraction.reached(specs, state),
    )
    @test !isempty(traj)
    @test traj[1] == initial_state
    @test length(traj) == length(ctrls) + 1
    @test all(x -> length(x) == 3, traj)
    for state in traj
        @test AB.TimedHybridAbstraction.reached(safety_specs, state)
    end
    println("Safety test completed: trajectory remained safe for $(length(traj)) steps")

    # Test minimal trajectory
    traj, ctrls = AB.TimedHybridAbstraction.get_closed_loop_trajectory(
        param_discretization,
        hs,
        safety_specs,
        concrete_controller,
        initial_state,
        1,
        stopping = (specs, state) -> !AB.TimedHybridAbstraction.reached(specs, state),
    )
    @test !isempty(traj)
    @test traj[1] == initial_state

    # Test safety checking function
    make_aug_state(xval, tval, kval) = ([xval], tval, kval)
    @test AB.TimedHybridAbstraction.reached(safety_specs, make_aug_state(0.0, 1.0, 1))
    @test AB.TimedHybridAbstraction.reached(safety_specs, make_aug_state(2.0, 3.0, 2))
    @test !AB.TimedHybridAbstraction.reached(safety_specs, make_aug_state(4.0, 1.0, 1))
    @test AB.TimedHybridAbstraction.reached(safety_specs, make_aug_state(-0.9, 2.0, 2))

    # Test solve shortcut
    controller = AB.TimedHybridAbstraction.solve_timed_hybrid_problem(
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
        safety_specs,
    )
    for state in traj[1:(end - 1)]
        @test controller.f(state) == concrete_controller.f(state)
    end
end

end # module
