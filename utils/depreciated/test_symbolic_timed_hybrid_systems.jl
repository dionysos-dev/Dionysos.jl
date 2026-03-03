module TestMain

using Test
using StaticArrays, Plots, MathematicalSystems, HybridSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

using MathOptInterface
const MOI = MathOptInterface

sleep(0.1)
println("Started test")

@testset "SymbolicTimedHybridSystems - simple, Time taken into account" begin
    #
    # Test 1: 2-mode hybrid system (1D state + 1D time)
    # - Mode 1: dx/dt = 0.5*x + u, x ∈ [-1,1], u ∈ [-0.5,0.5], time ∈ [0,2]
    # - Mode 2: dx/dt = 0.8*x + u, x ∈ [-2,2], u ∈ [-0.3,0.3], time ∈ [3,5]
    # - 1 transition (mode 1 → mode 2) with a reset map that forces x = 0.0, t = 3.0
    # - Reset map: FixedPointResetMap (ignores input state, always returns the target)
    #

    function mode1_dynamics(x, u)
        return [0.5 * x[1] + u[1]]
    end
    X1 = UT.HyperRectangle([-1.0], [1.0])
    U1 = UT.HyperRectangle([-0.5], [0.5])
    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode1_dynamics,
        1,
        1,
        X1,
        U1,
    )

    function mode2_dynamics(x, u)
        return [0.8 * x[1] + u[1]]
    end

    X2 = UT.HyperRectangle([-2.0], [2.0])
    U2 = UT.HyperRectangle([-0.3], [0.3])
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode2_dynamics,
        1,
        1,
        X2,
        U2,
    )

    time_mode1 = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([0.0], [2.0]),
    )

    time_mode2 = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([3.0], [5.0]),
    )

    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end

    function MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector)
        return reset.target
    end

    function MathematicalSystems.stateset(reset::FixedPointResetMap)
        return reset.domain
    end

    guard_1 = UT.HyperRectangle([0.0, 1.0], [1.0, 2.0])
    reset_map = FixedPointResetMap(guard_1, [0.0, 3.0])

    automaton = HybridSystems.GraphAutomaton(2)

    HybridSystems.add_transition!(automaton, 1, 2, 1)

    modes_systems = [
        SY.VectorContinuousSystem([mode1_system, time_mode1]),
        SY.VectorContinuousSystem([mode2_system, time_mode2]),
    ]
    reset_maps = [reset_map]

    switchings = [HybridSystems.AutonomousSwitching()]

    hs = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    state_grid_1 = DO.GridFree(SVector(0.0), SVector(0.5))
    state_grid_2 = DO.GridFree(SVector(0.0), SVector(0.5))
    input_grid_1 = DO.GridFree(SVector(0.0), SVector(0.5))
    input_grid_2 = DO.GridFree(SVector(0.0), SVector(0.3))

    optimizer_factory_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]
    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
    ]
    println(typeof(optimizer_factory_list[1]))
    # Using the new optimized function name
    hybrid_model = SY.SymbolicTimedHybridSystems.build_timed_hybrid_symbolic_model(
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
    )

    # ===================== TESTS =====================

    gim = hybrid_model.input_mapping
    @test gim.switching_inputs == 1
    @test length(gim.continuous_to_global) == 6
    @test length(gim.switching_to_global) == 1

    @test SY.SymbolicTimedHybridSystems.get_global_input_id(gim, 1, 1) == 1
    @test SY.SymbolicTimedHybridSystems.get_global_input_id(gim, 2, 1) == 4
    @test SY.SymbolicTimedHybridSystems.get_switching_global_id(gim, 1) == 7
    @test SY.SymbolicTimedHybridSystems.is_continuous_input(gim, 1) == true
    @test SY.SymbolicTimedHybridSystems.is_switching_input(
        gim,
        gim.switching_range.start,
    ) == true

    typ, info = SY.SymbolicTimedHybridSystems.get_local_input_info(gim, 1)
    @test typ == :continuous
    typ2, info2 =
        SY.SymbolicTimedHybridSystems.get_local_input_info(gim, gim.switching_range.start)
    @test typ2 == :switching

    for sm in hybrid_model.mode_abstractions
        @test SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_state(sm) > 0
        @test SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_input(sm) == 3
    end
    for tm in hybrid_model.time_abstractions
        @test length(tm.tsteps) > 0
        @test tm.tsteps[1] <= tm.tsteps[end]
    end

    autom = hybrid_model.symbolic_automaton
    for i in 1:length(hybrid_model.state_index_to_augmented)
        @test haskey(
            hybrid_model.augmented_to_state_index,
            hybrid_model.state_index_to_augmented[i],
        )
    end

    sm = hybrid_model.mode_abstractions[1]
    tm = hybrid_model.time_abstractions[1]
    x0 = [0.0]
    idx = SY.SymbolicTimedHybridSystems.find_symbolic_state(sm, x0)
    @test idx > 0
    guard = guard_1
    sp = SY.SymbolicTimedHybridSystems.extract_spatial_part(guard)
    tp = SY.SymbolicTimedHybridSystems.extract_temporal_part(guard)
    @test sp.lb == guard.lb[1:(end - 1)]
    @test tp == [guard.lb[end], guard.ub[end]]
    indices = SY.SymbolicTimedHybridSystems.get_time_indices_from_interval(tm, tp)
    @test all(i >= 1 && i <= length(tm.tsteps) for i in indices)
    idx_bad = SY.SymbolicTimedHybridSystems.find_symbolic_state(sm, [100.0])
    @test idx_bad == 0

    @test SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model) == 48

    @test SY.SymbolicTimedHybridSystems.get_n_input(hybrid_model) == 7

    @test collect(SY.SymbolicTimedHybridSystems.enum_states(hybrid_model)) ==
          collect(1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model))

    @test collect(SY.SymbolicTimedHybridSystems.enum_inputs(hybrid_model, 1)) ==
          collect(1:3)
    @test collect(SY.SymbolicTimedHybridSystems.enum_inputs(hybrid_model, 2)) ==
          collect(4:6) .- 3

    for state in 1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model)
        aug = hybrid_model.state_index_to_augmented[state]
        concrete = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_model, state)
        idx = SY.SymbolicTimedHybridSystems.get_abstract_state(hybrid_model, concrete)
        @test idx == state
    end

    Xs = [UT.HyperRectangle([-1.0], [1.0]), UT.HyperRectangle([-2.0], [2.0])]
    Ts = [UT.HyperRectangle([0.0], [2.0]), UT.HyperRectangle([3.0], [5.0])]
    Ns = [1]

    # Appel de la fonction
    indices = SY.SymbolicTimedHybridSystems.get_states_from_set(hybrid_model, Xs, Ts, Ns)

    all_mode1 = [
        i for i in 1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model) if
        hybrid_model.state_index_to_augmented[i][3] == 1
    ]
    @test sort(indices) == sort(all_mode1)

    Xs2 = [UT.HyperRectangle([0.0], [0.5]), UT.HyperRectangle([-2.0], [2.0])]
    Ts2 = [UT.HyperRectangle([0.0], [1.0]), UT.HyperRectangle([3.0], [5.0])]
    indices2 = SY.SymbolicTimedHybridSystems.get_states_from_set(hybrid_model, Xs2, Ts2, Ns)
    for idx in indices2
        aug = hybrid_model.state_index_to_augmented[idx]
        @test aug[3] == 1
        tval = hybrid_model.time_abstractions[1].tsteps[aug[2]]
        @test 0.0 <= tval <= 1.0
        xval = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_model, idx)[1][1]
        @test 0.0 <= xval <= 0.5
    end
end

@testset "SymbolicTimedHybridSystems - simple, Time not taken into account" begin
    #
    # Test 2: 2-mode hybrid system (1D state + 1D time, frozen time)
    # - Mode 1: dx/dt = 0.5*x + u, x ∈ [-1,1], u ∈ [-0.5,0.5], time ∈ [0,2]
    # - Mode 2: dx/dt = 0.8*x + u, x ∈ [-2,2], u ∈ [-0.3,0.3], time ∈ [3,5]
    # - 1 transition (mode 1 → mode 2) with a reset map that forces x = 0.0, t = 0.0
    # - Reset map: FixedPointResetMap (ignores input state, always returns the target)
    #
    function mode1_dynamics(x, u)
        return [0.5 * x[1] + u[1]]
    end
    X1 = UT.HyperRectangle([-1.0], [1.0])
    U1 = UT.HyperRectangle([-0.5], [0.5])
    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode1_dynamics,
        1,
        1,
        X1,
        U1,
    )

    function mode2_dynamics(x, u)
        return [0.8 * x[1] + u[1]]
    end
    X2 = UT.HyperRectangle([-2.0], [2.0])
    U2 = UT.HyperRectangle([-0.3], [0.3])
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode2_dynamics,
        1,
        1,
        X2,
        U2,
    )

    temps_mode1 = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [0.0;;],
        UT.HyperRectangle([0.0], [2.0]),
    )
    temps_mode2 = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [0.0;;],
        UT.HyperRectangle([3.0], [5.0]),
    )

    struct FixedPointResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle
        target::Vector{Float64}
    end

    function MathematicalSystems.apply(reset::FixedPointResetMap, state::AbstractVector)
        return reset.target
    end

    function MathematicalSystems.stateset(reset::FixedPointResetMap)
        return reset.domain
    end

    guard_1 = UT.HyperRectangle([0.0, -1.0], [1.0, 1.0])  # Guard: x ∈ [0,1], t ∈ [-1,1] by default
    reset_map = FixedPointResetMap(guard_1, [0.0, 0.0])

    automaton = HybridSystems.GraphAutomaton(2)

    HybridSystems.add_transition!(automaton, 1, 2, 1)

    modes_systems = [
        SY.VectorContinuousSystem([mode1_system, temps_mode1]),
        SY.VectorContinuousSystem([mode2_system, temps_mode2]),
    ]

    reset_maps = [reset_map]
    switchings = [HybridSystems.AutonomousSwitching()]

    hs = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Configuration des grilles et optimizers pour chaque mode
    state_grid_1 = DO.GridFree(SVector(0.0), SVector(0.5))
    state_grid_2 = DO.GridFree(SVector(0.0), SVector(0.5))
    input_grid_1 = DO.GridFree(SVector(0.0), SVector(0.5))
    input_grid_2 = DO.GridFree(SVector(0.0), SVector(0.3))

    optimizer_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]
    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.8),
        ),
    ]

    hybrid_model = SY.SymbolicTimedHybridSystems.build_timed_hybrid_symbolic_model(
        hs,
        optimizer_list,
        optimizer_kwargs_dict,
    )

    # ===================== TESTS =====================

    for aug_state in hybrid_model.state_index_to_augmented
        # aug_state = (state_symbol, time_symbol, mode_id)
        @test aug_state[2] == 1
    end

    for tm in hybrid_model.time_abstractions
        @test length(tm.tsteps) == 1
        @test tm.is_time_active == false
    end

    for tm in hybrid_model.time_abstractions
        tmin, tmax = tm.tsteps[1], tm.tsteps[1]
        @test SY.SymbolicTimedHybridSystems.get_time_indices_from_interval(
            tm,
            [tmin, tmax],
        ) == [1]
    end

    @test SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model) ==
          length(hybrid_model.state_index_to_augmented)
    @test collect(SY.SymbolicTimedHybridSystems.enum_states(hybrid_model)) ==
          collect(1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model))

    @test SY.SymbolicTimedHybridSystems.get_n_input(hybrid_model) ==
          hybrid_model.input_mapping.total_inputs

    for k in 1:length(hybrid_model.mode_abstractions)
        expected = collect(
            1:SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_input(
                hybrid_model.mode_abstractions[k],
            ),
        )
        actual = collect(SY.SymbolicTimedHybridSystems.enum_inputs(hybrid_model, k))
        @test actual == expected
    end

    for state in 1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model)
        aug = hybrid_model.state_index_to_augmented[state]
        concrete = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_model, state)
        idx = SY.SymbolicTimedHybridSystems.get_abstract_state(hybrid_model, concrete)
        @test idx == state
    end

    gim = hybrid_model.input_mapping
    for k in 1:length(hybrid_model.mode_abstractions)
        n_inputs = SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_input(
            hybrid_model.mode_abstractions[k],
        )
        for local_input_id in 1:n_inputs
            global_input_id =
                SY.SymbolicTimedHybridSystems.get_global_input_id(gim, k, local_input_id)

            u_concrete = SY.SymbolicTimedHybridSystems.get_concrete_input(
                hybrid_model,
                global_input_id,
                k,
            )
            global_input_id2 = SY.SymbolicTimedHybridSystems.get_abstract_input(
                hybrid_model,
                u_concrete,
                k,
            )
            @test global_input_id == global_input_id2
        end
    end

    for tr in 1:gim.switching_inputs
        gid = SY.SymbolicTimedHybridSystems.get_switching_global_id(gim, tr)
        @test isnothing(
            SY.SymbolicTimedHybridSystems.get_concrete_input(hybrid_model, gid, 1),
        )
    end

    if hasproperty(gim, :switch_labels)
        @test length(gim.switch_labels) == gim.switching_inputs
        for label in gim.switch_labels
            @test occursin("SWITCH", label)
            @test occursin("->", label)
        end
    end
end

@testset "SymbolicTimedHybridSystems - Complex Multi-Mode System" begin
    # Test case: 3-mode hybrid system with complex dynamics
    # Mode 1: Linear system
    # Mode 2: Affine system  
    # Mode 3: Black box system
    # Multiple transitions with different guards and resets

    # Mode 1: Simple linear dynamics
    function mode1_dynamics(x, u)
        return [0.3 * x[1] + 0.1 * u[1]]
    end
    X1 = UT.HyperRectangle([-2.0], [2.0])
    U1 = UT.HyperRectangle([-1.0], [1.0])
    mode1_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode1_dynamics,
        1,
        1,
        X1,
        U1,
    )

    # Mode 2: More complex nonlinear dynamics
    function mode2_dynamics(x, u)
        return [0.5 * x[1]^2 * sign(x[1]) + u[1]]
    end
    X2 = UT.HyperRectangle([-3.0], [3.0])
    U2 = UT.HyperRectangle([-2.0], [2.0])
    mode2_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode2_dynamics,
        1,
        1,
        X2,
        U2,
    )

    # Mode 3: Oscillatory dynamics
    function mode3_dynamics(x, u)
        return [-0.2 * x[1] + 0.1 * sin(x[1]) + u[1]]
    end
    X3 = UT.HyperRectangle([-1.5], [1.5])
    U3 = UT.HyperRectangle([-0.5], [0.5])
    mode3_system = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        mode3_dynamics,
        1,
        1,
        X3,
        U3,
    )

    # Time dynamics for each mode (different time evolutions)
    time_mode1 = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([0.0], [5.0]),
    )
    time_mode2 = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([2.0], [8.0]),
    )
    time_mode3 = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([1.0], [6.0]),
    )

    # Define reset maps
    struct LinearResetMap <: MathematicalSystems.AbstractMap
        A::Matrix{Float64}
        b::Vector{Float64}
        domain::UT.HyperRectangle
    end

    function MathematicalSystems.apply(reset::LinearResetMap, state::AbstractVector)
        return reset.A * state + reset.b
    end

    function MathematicalSystems.stateset(reset::LinearResetMap)
        return reset.domain
    end

    # Create reset maps for transitions
    # Transition 1→2: Scale and shift
    guard_1_2 = UT.HyperRectangle([1.0, 3.0], [2.0, 5.0])
    reset_1_2 = LinearResetMap([0.8 0.0; 0.0 1.0], [0.5, 2.5], guard_1_2)

    # Transition 2→3: More complex transformation
    guard_2_3 = UT.HyperRectangle([2.0, 5.0], [3.0, 8.0])
    reset_2_3 = LinearResetMap([0.5 0.0; 0.0 0.8], [-0.2, 2.0], guard_2_3)

    # Transition 3→1: Reset to center
    guard_3_1 = UT.HyperRectangle([-1.0, 4.0], [1.0, 6.0])
    reset_3_1 = LinearResetMap([0.1 0.0; 0.0 0.2], [0.0, 1.0], guard_3_1)

    # Build automaton with 3 modes and 3 transitions
    automaton = HybridSystems.GraphAutomaton(3)
    HybridSystems.add_transition!(automaton, 1, 2, 1)  # 1→2
    HybridSystems.add_transition!(automaton, 2, 3, 2)  # 2→3 
    HybridSystems.add_transition!(automaton, 3, 1, 3)  # 3→1

    modes_systems = [
        SY.VectorContinuousSystem([mode1_system, time_mode1]),
        SY.VectorContinuousSystem([mode2_system, time_mode2]),
        SY.VectorContinuousSystem([mode3_system, time_mode3]),
    ]
    reset_maps = [reset_1_2, reset_2_3, reset_3_1]
    switchings = [
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
        HybridSystems.AutonomousSwitching(),
    ]

    hs = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Configure different grid sizes for each mode
    state_grid_1 = DO.GridFree(SVector(0.0), SVector(0.4))  # Fine grid
    state_grid_2 = DO.GridFree(SVector(0.0), SVector(0.6))  # Medium grid
    state_grid_3 = DO.GridFree(SVector(0.0), SVector(0.3))  # Very fine grid

    input_grid_1 = DO.GridFree(SVector(0.0), SVector(0.5))
    input_grid_2 = DO.GridFree(SVector(0.0), SVector(0.8))
    input_grid_3 = DO.GridFree(SVector(0.0), SVector(0.25))

    optimizer_factory_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_1,
            "input_grid" => input_grid_1,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.4),
        ),
        Dict(
            "state_grid" => state_grid_2,
            "input_grid" => input_grid_2,
            "time_step" => 0.6,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(1.0),
        ),
        Dict(
            "state_grid" => state_grid_3,
            "input_grid" => input_grid_3,
            "time_step" => 0.4,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.3),
        ),
    ]

    # Build the hybrid model
    hybrid_model = SY.SymbolicTimedHybridSystems.build_timed_hybrid_symbolic_model(
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
    )

    # ===================== COMPREHENSIVE TESTS =====================

    # Test basic structure
    @test length(hybrid_model.mode_abstractions) == 3
    @test length(hybrid_model.time_abstractions) == 3
    @test hybrid_model.input_mapping.switching_inputs == 3

    # Test that each mode has the expected number of inputs
    for (i, mode_abs) in enumerate(hybrid_model.mode_abstractions)
        n_inputs = Dionysos.Symbolic.get_n_input(mode_abs)
        @test n_inputs > 0  # Should have at least some inputs
    end

    # Test time abstractions are correct
    for (i, time_abs) in enumerate(hybrid_model.time_abstractions)
        @test length(time_abs.tsteps) > 1  # Should have multiple time steps
        @test time_abs.is_time_active == true
    end

    # Test state mappings consistency
    n_states = SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model)
    @test n_states == length(hybrid_model.state_index_to_augmented)
    @test n_states == length(hybrid_model.augmented_to_state_index)

    # Test that all states can be mapped back and forth
    for state_idx in 1:min(n_states, 100)  # Test first 100 states
        aug_state = hybrid_model.state_index_to_augmented[state_idx]
        @test haskey(hybrid_model.augmented_to_state_index, aug_state)
        @test hybrid_model.augmented_to_state_index[aug_state] == state_idx
    end

    # Test input mappings
    gim = hybrid_model.input_mapping
    @test gim.total_inputs == gim.continuous_inputs + gim.switching_inputs
    @test gim.continuous_inputs > 0
    @test gim.switching_inputs == 3  # One for each transition

    # Test switching input labels
    @test length(gim.switch_labels) == 3
    for label in gim.switch_labels
        @test occursin("SWITCH", label)
        @test occursin("->", label)
    end

    # Test concrete/abstract state conversions for a sample of states
    sample_states = min(20, n_states)
    for state_idx in 1:sample_states
        concrete = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_model, state_idx)
        @test length(concrete) == 3  # (continuous_state, time, mode_id)

        # Verify mode_id is valid
        mode_id = concrete[3]
        @test 1 <= mode_id <= 3
    end

    println("✓ Complex multi-mode system test completed with $(n_states) states")
end

@testset "SymbolicTimedHybridSystems - Performance and Memory Tests" begin
    # Test performance characteristics of the optimized functions

    # Create a moderately complex system for performance testing
    function dynamics(x, u)
        return [0.5 * x[1] + u[1]]
    end

    X = UT.HyperRectangle([-2.0], [2.0])
    U = UT.HyperRectangle([-1.0], [1.0])
    system =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(dynamics, 1, 1, X, U)

    time_system = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([0.0], [3.0]),
    )

    struct SimpleResetMap <: MathematicalSystems.AbstractMap
        target::Vector{Float64}
        domain::UT.HyperRectangle
    end

    function MathematicalSystems.apply(reset::SimpleResetMap, state::AbstractVector)
        return reset.target
    end

    function MathematicalSystems.stateset(reset::SimpleResetMap)
        return reset.domain
    end

    guard = UT.HyperRectangle([1.0, 2.0], [2.0, 3.0])
    reset_map = SimpleResetMap([0.0, 0.5], guard)

    automaton = HybridSystems.GraphAutomaton(2)
    HybridSystems.add_transition!(automaton, 1, 2, 1)

    modes_systems = [
        SY.VectorContinuousSystem([system, time_system]),
        SY.VectorContinuousSystem([system, time_system]),
    ]
    reset_maps = [reset_map]
    switchings = [HybridSystems.AutonomousSwitching()]

    hs = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Use fine grids to create a larger system for performance testing
    state_grid = DO.GridFree(SVector(0.0), SVector(0.2))  # Fine grid
    input_grid = DO.GridFree(SVector(0.0), SVector(0.2))

    optimizer_factory_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]

    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid,
            "input_grid" => input_grid,
            "time_step" => 0.3,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
        Dict(
            "state_grid" => state_grid,
            "input_grid" => input_grid,
            "time_step" => 0.3,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{1, 1}(0.5),
        ),
    ]

    # Measure construction time
    construction_time = @elapsed begin
        hybrid_model = SY.SymbolicTimedHybridSystems.build_timed_hybrid_symbolic_model(
            hs,
            optimizer_factory_list,
            optimizer_kwargs_dict,
        )
    end

    n_states = SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model)
    n_inputs = SY.SymbolicTimedHybridSystems.get_n_input(hybrid_model)

    println(
        "✓ Performance test: Built system with $(n_states) states, $(n_inputs) inputs in $(round(construction_time, digits=3))s",
    )

    # Test that construction time is reasonable (less than 10 seconds for this size)
    @test construction_time < 10.0

    # Test basic functionality works correctly even with larger systems
    @test n_states > 0
    @test n_inputs > 0

    # Test a sample of state conversions for correctness
    sample_size = min(50, n_states)
    conversion_time = @elapsed begin
        for state_idx in 1:sample_size
            concrete = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_model, state_idx)
            @test length(concrete) == 3  # (continuous_state, time, mode_id)
        end
    end

    # Test that state conversion is efficient
    avg_conversion_time = conversion_time / sample_size
    @test avg_conversion_time < 0.01  # Should be very fast per conversion

    println(
        "✓ State conversion performance: $(round(avg_conversion_time * 1000, digits=3))ms per conversion",
    )
end

@testset "SymbolicTimedHybridSystems - 3 modes, transitions complètes, reset maps variées" begin
    #
    # Test 4: 3-mode hybrid system (2D state + 1D time)
    # - Mode A: dx/dt = A_A*x + B_A*u, x ∈ [-2,2]^2, u ∈ [-1,1]^2, time ∈ [0,3]
    # - Mode B: dx/dt = A_B*x + B_B*u, x ∈ [-2,2]^2, u ∈ [-2,1]^2, time ∈ [1.5,2.5]
    # - Mode C: dx/dt = A_C*x + B_C*u, x ∈ [-2,2]^2, u ∈ [-1,1]^2, time ∈ [2,3]
    # - 6 transitions (all directions) with various reset maps:
    #   * offset (translation), perm (component permutation), tmin (lower bound on t)
    #   * each reset map bounds x_new by the upper bound of the target mode
    # - Reset map: DenseResetMap (permutation, offset, tmin, x_max)
    #

    A_A = [0.5 0.0; 0.0 0.6]
    B_A = [1.0 0.0; 0.0 1.0]
    A_B = [0.4 0.1; 0.1 0.7]
    B_B = [1.0 0.0; 0.0 1.0]
    A_C = [0.8 0.1; 0.1 0.9]
    B_C = [1.0 0.0; 0.2 1.0]

    X_A = UT.HyperRectangle([-2.0, -2.0], [2.0, 2.0])
    U_A = UT.HyperRectangle([-1.0, -1.0], [1.0, 1.0])
    X_B = UT.HyperRectangle([-2.0, -2.0], [2.0, 2.0])
    U_B = UT.HyperRectangle([-2.0, -2.0], [1.0, 1.0])
    X_C = UT.HyperRectangle([-2.0, -2.0], [2.0, 2.0])
    U_C = UT.HyperRectangle([-1.0, -1.0], [1.0, 1.0])

    function dynA(x, u)
        return A_A * x + B_A * u
    end
    function dynB(x, u)
        return A_B * x + B_B * u
    end
    function dynC(x, u)
        return A_C * x + B_C * u
    end
    modeA =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(dynA, 2, 2, X_A, U_A)
    modeB =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(dynB, 2, 2, X_B, U_B)
    modeC =
        MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(dynC, 2, 2, X_C, U_C)

    temps_A = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([0.0], [3.0]),
    )
    temps_B = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([1.5], [2.5]),
    )
    temps_C = MathematicalSystems.ConstrainedLinearContinuousSystem(
        [1.0;;],
        UT.HyperRectangle([2.0], [3.0]),
    )

    vcA = SY.VectorContinuousSystem([modeA, temps_A])
    vcB = SY.VectorContinuousSystem([modeB, temps_B])
    vcC = SY.VectorContinuousSystem([modeC, temps_C])

    automaton = HybridSystems.GraphAutomaton(3)
    HybridSystems.add_transition!(automaton, 1, 2, 1) # A→B
    HybridSystems.add_transition!(automaton, 2, 1, 2) # B→A
    HybridSystems.add_transition!(automaton, 1, 3, 3) # A→C
    HybridSystems.add_transition!(automaton, 3, 1, 4) # C→A
    HybridSystems.add_transition!(automaton, 2, 3, 5) # B→C
    HybridSystems.add_transition!(automaton, 3, 2, 6) # C→B

    # === Various reset maps ===
    # A reset map defines how the state changes during a transition between two modes.
    # Here, DenseResetMap allows modeling complex resets:
    #  - domain : domain of validity of the reset (guard associated with the transition)
    #  - offset : translation applied to the state part (x)
    #  - perm   : permutation of the state components (e.g., [2,1] swaps x₁ and x₂)
    #  - tmin   : minimal value imposed on the time component after reset
    struct DenseResetMap <: MathematicalSystems.AbstractMap
        domain::UT.HyperRectangle  # Domain on which the reset is valid (guard)
        offset::Vector{Float64}    # Offset applied to the state (x)
        perm::Vector{Int}          # Permutation of the state indices (x)
        tmin::Float64              # Minimal time after reset
        x_max::Vector{Float64}     # Maximum allowed value for each component of x after reset (upper bound of the target mode)
    end

    # Reset map application function:
    # - Takes the current state (x₁, x₂, t)
    # - Applies the permutation perm to x
    # - Adds the offset to the permuted x
    # - Updates t by enforcing t >= tmin
    # - Returns the new state [x_new; t_new]
    function MathematicalSystems.apply(reset::DenseResetMap, state::AbstractVector)
        x = state[1:2]  # State part (assumed dimension 2)
        t = max(state[3], reset.tmin)  # Time bounded below by tmin
        x_perm = x[reset.perm]         # Permute components
        # Apply offset then bound by x_max (component-wise)
        x_new = min.(x_perm .+ reset.offset, reset.x_max)
        return vcat(x_new, t)          # Concatenate new state
    end

    function MathematicalSystems.stateset(reset::DenseResetMap)
        return reset.domain
    end

    guards = [
        UT.HyperRectangle([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),   # A→B : x₁∈[1,2], x₂∈[0,1], t∈[0,1]
        UT.HyperRectangle([0.0, 1.0, 1.5], [1.0, 2.0, 2.5]),   # B→A : x₁∈[0,1], x₂∈[1,2], t∈[1.5,2.5]
        UT.HyperRectangle([1.0, -1.0, 0.0], [2.0, 0.0, 1.0]),  # A→C : x₁∈[1,2], x₂∈[-1,0], t∈[0,1]
        UT.HyperRectangle([-1.0, 1.0, 2.0], [0.0, 2.0, 3.0]),  # C→A : x₁∈[-1,0], x₂∈[1,2], t∈[2,3]
        UT.HyperRectangle([1.0, 1.0, 1.5], [2.0, 2.0, 2.5]),   # B→C : x₁∈[1,2], x₂∈[1,2], t∈[1.5,2.5]
        UT.HyperRectangle([-1.0, -1.0, 2.0], [0.0, 0.0, 3.0]),  # C→B : x₁∈[-1,0], x₂∈[-1,0], t∈[2,3]
    ]

    offsets = [
        [0.1, 0.0],   # A→B : translation
        [0.0, -0.1],  # B→A : translation négative
        [0.0, 0.0],   # A→C : identité
        [0.0, 0.2],  # C→A : translation mixte
        [0.0, 0.0],   # B→C : identité
        [0.3, 0.0],    # C→B : translation positive
    ]
    perms = [
        [1, 2],  # A→B : identity
        [2, 1],  # B→A : permutation
        [1, 2],  # A→C : identity
        [2, 1],  # C→A : permutation
        [1, 2],  # B→C : identity
        [2, 1],  # C→B : permutation
    ]
    tmins = [1.5, 0.0, 2.0, 0.0, 2.0, 1.5]

    x_maxs = [X_B.ub, X_A.ub, X_C.ub, X_A.ub, X_C.ub, X_B.ub]  # transitions : 1→2, 2→1, 1→3, 3→1, 2→3, 3→2
    reset_maps = [
        DenseResetMap(guards[1], offsets[1], perms[1], tmins[1], x_maxs[1]),
        DenseResetMap(guards[2], offsets[2], perms[2], tmins[2], x_maxs[2]),
        DenseResetMap(guards[3], offsets[3], perms[3], tmins[3], x_maxs[3]),
        DenseResetMap(guards[4], offsets[4], perms[4], tmins[4], x_maxs[4]),
        DenseResetMap(guards[5], offsets[5], perms[5], tmins[5], x_maxs[5]),
        DenseResetMap(guards[6], offsets[6], perms[6], tmins[6], x_maxs[6]),
    ]
    switchings = [HybridSystems.AutonomousSwitching() for _ in 1:6]

    modes_systems = [vcA, vcB, vcC]

    hs = HybridSystems.HybridSystem(automaton, modes_systems, reset_maps, switchings)

    # Configuration des grilles et optimizers pour chaque mode (3 modes)
    state_grid_A = DO.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5))
    state_grid_B = DO.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5))
    state_grid_C = DO.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5))
    input_grid_A = DO.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5))
    input_grid_B = DO.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5))
    input_grid_C = DO.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5))

    optimizer_factory_list = [
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
        () -> MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer),
    ]
    optimizer_kwargs_dict = [
        Dict(
            "state_grid" => state_grid_A,
            "input_grid" => input_grid_A,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{2, 2}(abs.(A_A)),
        ),
        Dict(
            "state_grid" => state_grid_B,
            "input_grid" => input_grid_B,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{2, 2}(abs.(A_B)),
        ),
        Dict(
            "state_grid" => state_grid_C,
            "input_grid" => input_grid_C,
            "time_step" => 0.5,
            "approx_mode" => Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            "jacobian_bound" => u -> SMatrix{2, 2}(abs.(A_C)),
        ),
    ]

    # Using the new optimized API
    hybrid_model = SY.SymbolicTimedHybridSystems.build_timed_hybrid_symbolic_model(
        hs,
        optimizer_factory_list,
        optimizer_kwargs_dict,
    )

    # ===================== TESTS =====================

    # Test basic structure with new API
    @test length(hybrid_model.mode_abstractions) == 3
    @test length(hybrid_model.time_abstractions) == 3
    @test length(hybrid_model.state_index_to_augmented) ==
          length(hybrid_model.augmented_to_state_index)

    gim = hybrid_model.input_mapping
    @test gim.switching_inputs == 6
    @test length(gim.continuous_to_global) == 25 + 49 + 25  # mode A: 25, mode B: 49, mode C: 25
    @test length(gim.switching_to_global) == 6

    # Test continuous input mapping
    for mode in 1:3
        for local_input in
            1:min(
            5,
            SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_input(
                hybrid_model.mode_abstractions[mode],
            ),
        )  # Test with first few inputs

            gid = SY.SymbolicTimedHybridSystems.get_global_input_id(gim, mode, local_input)
            @test gid > 0
            typ, info = SY.SymbolicTimedHybridSystems.get_local_input_info(gim, gid)
            @test typ == :continuous
            @test info[1] == mode
        end
    end

    # Test switching input mapping
    for tr in 1:6
        gid = SY.SymbolicTimedHybridSystems.get_switching_global_id(gim, tr)
        @test gid > 0
        typ, info = SY.SymbolicTimedHybridSystems.get_local_input_info(gim, gid)
        @test typ == :switching
        @test info == tr
    end

    # Test mode abstractions
    mode_n_input = [25, 49, 25]  # Expected based on uniform grids
    for (idxsm, sm) in enumerate(hybrid_model.mode_abstractions)
        @test SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_state(sm) > 0
        @test SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_input(sm) ==
              mode_n_input[idxsm]
    end

    # Test time abstractions
    for tm in hybrid_model.time_abstractions
        @test length(tm.tsteps) > 0
        @test tm.tsteps[1] <= tm.tsteps[end]
    end

    # Test state mapping consistency
    for i in 1:length(hybrid_model.state_index_to_augmented)
        @test haskey(
            hybrid_model.augmented_to_state_index,
            hybrid_model.state_index_to_augmented[i],
        )
    end

    # Test symbolic state finding and guard processing
    for mode in 1:3
        sm = hybrid_model.mode_abstractions[mode]
        tm = hybrid_model.time_abstractions[mode]
        x0 = [0.0, 0.0]
        idx = SY.SymbolicTimedHybridSystems.find_symbolic_state(sm, x0)
        @test idx > 0

        idx_bad = SY.SymbolicTimedHybridSystems.find_symbolic_state(sm, [100.0, 100.0])
        @test idx_bad == 0

        guard = guards[mode]
        sp = SY.SymbolicTimedHybridSystems.extract_spatial_part(guard)
        tp = SY.SymbolicTimedHybridSystems.extract_temporal_part(guard)
        @test sp.lb == guard.lb[1:(end - 1)]
        @test tp == [guard.lb[end], guard.ub[end]]
        indices = SY.SymbolicTimedHybridSystems.get_time_indices_from_interval(tm, tp)
        @test all(i >= 1 && i <= length(tm.tsteps) for i in indices)
    end

    # Test invalid state handling
    for mode in 1:3
        sm = hybrid_model.mode_abstractions[mode]
        @test SY.SymbolicTimedHybridSystems.find_symbolic_state(sm, [999.0, -999.0]) == 0
    end

    # Test global input info consistency
    for gid in 1:gim.total_inputs
        typ, info = SY.SymbolicTimedHybridSystems.get_local_input_info(gim, gid)
        if typ == :continuous
            @test 1 <= info[1] <= 3  # mode_id should be 1, 2, or 3
            @test 1 <= info[2] <= 49  # local_input_id should be within valid range (max 49 for mode B)
        elseif typ == :switching
            @test 1 <= info <= 6
        else
            @test typ == :invalid
        end
    end

    # Test augmented state validity
    for i in 1:length(hybrid_model.state_index_to_augmented)
        aug = hybrid_model.state_index_to_augmented[i]
        mode = aug[3]
        t_idx = aug[2]
        @test 1 <= mode <= 3
        @test 1 <= t_idx <= length(hybrid_model.time_abstractions[mode].tsteps)
    end

    # Test accessors with new API
    @test SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model) ==
          length(hybrid_model.state_index_to_augmented)
    @test SY.SymbolicTimedHybridSystems.get_n_input(hybrid_model) ==
          hybrid_model.input_mapping.total_inputs

    @test collect(SY.SymbolicTimedHybridSystems.enum_states(hybrid_model)) ==
          collect(1:SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model))

    # Test input enumeration per mode
    for k in 1:length(hybrid_model.mode_abstractions)
        expected = collect(
            1:SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_input(
                hybrid_model.mode_abstractions[k],
            ),
        )
        actual = collect(SY.SymbolicTimedHybridSystems.enum_inputs(hybrid_model, k))
        @test actual == expected
    end

    # Test state conversion consistency
    for state in 1:min(100, SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model))  # Test first 100 states
        aug = hybrid_model.state_index_to_augmented[state]
        concrete = SY.SymbolicTimedHybridSystems.get_concrete_state(hybrid_model, state)
        idx = SY.SymbolicTimedHybridSystems.get_abstract_state(hybrid_model, concrete)
        @test idx == state
    end

    # Test input conversion consistency
    gim = hybrid_model.input_mapping
    for k in 1:length(hybrid_model.mode_abstractions)
        n_inputs = SY.SymbolicTimedHybridSystems.Dionysos.Symbolic.get_n_input(
            hybrid_model.mode_abstractions[k],
        )
        for local_input_id in 1:min(5, n_inputs)  # Test first few inputs
            global_input_id =
                SY.SymbolicTimedHybridSystems.get_global_input_id(gim, k, local_input_id)
            u_concrete = SY.SymbolicTimedHybridSystems.get_concrete_input(
                hybrid_model,
                global_input_id,
                k,
            )
            global_input_id2 = SY.SymbolicTimedHybridSystems.get_abstract_input(
                hybrid_model,
                u_concrete,
                k,
            )
            @test global_input_id == global_input_id2
        end
    end

    # Test switching inputs return nothing for concrete input
    for tr in 1:gim.switching_inputs
        gid = SY.SymbolicTimedHybridSystems.get_switching_global_id(gim, tr)
        @test isnothing(
            SY.SymbolicTimedHybridSystems.get_concrete_input(hybrid_model, gid, 1),
        )
    end

    # Test switching input labels
    expected_labels = [
        "SWITCH 1 -> 2",
        "SWITCH 2 -> 1",
        "SWITCH 1 -> 3",
        "SWITCH 3 -> 1",
        "SWITCH 2 -> 3",
        "SWITCH 3 -> 2",
    ]
    @test length(gim.switch_labels) == length(expected_labels)
    for label in expected_labels
        @test label in gim.switch_labels
    end
    for label in gim.switch_labels
        @test occursin("SWITCH", label)
        @test occursin("->", label)
    end

    println(
        "✓ Complex 3-mode system with varied reset maps: $(SY.SymbolicTimedHybridSystems.get_n_state(hybrid_model)) states",
    )
end

println("End test")
end
