module TestMain

using Test
using StaticArrays, Plots, MathematicalSystems, HybridSystems
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

sleep(0.1)
println("Started test")

@testset "TimedHybridAutomata - simple, Time taken into account" begin
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

    growth_bounds = SVector(SMatrix{1, 1}(0.5), SMatrix{1, 1}(0.8))

    param_discretization = [(0.5, 0.5, 0.5), (0.5, 0.3, 0.5)]
    hybrid_model = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        growth_bounds,
        param_discretization,
    )

    # ===================== TESTS =====================

    gim = hybrid_model.global_input_map
    @test gim.switching_inputs == 1
    @test length(gim.continuous_to_global) == 6
    @test length(gim.switching_to_global) == 1

    @test SY.TimedHybridAutomata.get_global_input_id(gim, 1, 1) == 1
    @test SY.TimedHybridAutomata.get_global_input_id(gim, 2, 1) == 4
    @test SY.TimedHybridAutomata.get_switching_global_id(gim, 1) == 7
    @test SY.TimedHybridAutomata.is_continuous_input(gim, 1) == true
    @test SY.TimedHybridAutomata.is_switching_input(gim, gim.switching_range.start) == true

    typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, 1)
    @test typ == :continuous
    typ2, info2 =
        SY.TimedHybridAutomata.get_local_input_info(gim, gim.switching_range.start)
    @test typ2 == :switching

    for sm in hybrid_model.symmodels
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_state(sm) > 0
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(sm) == 3
    end
    for tm in hybrid_model.time_symbolic_models
        @test length(tm.tsteps) > 0
        @test tm.tsteps[1] <= tm.tsteps[end]
    end

    autom = hybrid_model.autom
    for i in 1:length(hybrid_model.int2aug_state)
        @test haskey(hybrid_model.aug_state2int, hybrid_model.int2aug_state[i])
    end

    sm = hybrid_model.symmodels[1]
    tm = hybrid_model.time_symbolic_models[1]
    x0 = [0.0]
    idx = SY.TimedHybridAutomata.find_symbolic_state(sm, x0)
    @test idx > 0
    t0 = tm.tsteps[1]
    tidx = SY.TimedHybridAutomata.find_time_index(tm, t0)
    @test tidx == 1
    guard = guard_1
    sp = SY.TimedHybridAutomata.extract_spatial_part(guard)
    tp = SY.TimedHybridAutomata.extract_temporal_part(guard)
    @test sp.lb == guard.lb[1:(end - 1)]
    @test tp == [guard.lb[end], guard.ub[end]]
    indices = SY.TimedHybridAutomata.get_time_indices_from_interval(tm, tp)
    @test all(i >= 1 && i <= length(tm.tsteps) for i in indices)
    idx_bad = SY.TimedHybridAutomata.find_symbolic_state(sm, [100.0])
    @test idx_bad == 0

    @test SY.TimedHybridAutomata.get_n_state(hybrid_model) == 48

    @test SY.TimedHybridAutomata.get_n_input(hybrid_model) == 7

    @test collect(SY.TimedHybridAutomata.enum_states(hybrid_model)) ==
          collect(1:SY.TimedHybridAutomata.get_n_state(hybrid_model))

    @test collect(SY.TimedHybridAutomata.enum_inputs(hybrid_model, 1)) == collect(1:3)
    @test collect(SY.TimedHybridAutomata.enum_inputs(hybrid_model, 2)) == collect(4:6) .- 3

    for state in 1:SY.TimedHybridAutomata.get_n_state(hybrid_model)
        aug = hybrid_model.int2aug_state[state]
        concrete = SY.TimedHybridAutomata.get_concrete_state(hybrid_model, state)
        idx = SY.TimedHybridAutomata.get_abstract_state(hybrid_model, concrete)
        @test idx == state
    end

    Xs = [UT.HyperRectangle([-1.0], [1.0]), UT.HyperRectangle([-2.0], [2.0])]
    Ts = [UT.HyperRectangle([0.0], [2.0]), UT.HyperRectangle([3.0], [5.0])]
    Ns = [1]

    # Appel de la fonction
    indices = SY.TimedHybridAutomata.get_states_from_set(hybrid_model, Xs, Ts, Ns)

    all_mode1 = [
        i for i in 1:SY.TimedHybridAutomata.get_n_state(hybrid_model) if
        hybrid_model.int2aug_state[i][3] == 1
    ]
    @test sort(indices) == sort(all_mode1)

    Xs2 = [UT.HyperRectangle([0.0], [0.5]), UT.HyperRectangle([-2.0], [2.0])]
    Ts2 = [UT.HyperRectangle([0.0], [1.0]), UT.HyperRectangle([3.0], [5.0])]
    indices2 = SY.TimedHybridAutomata.get_states_from_set(hybrid_model, Xs2, Ts2, Ns)
    for idx in indices2
        aug = hybrid_model.int2aug_state[idx]
        @test aug[3] == 1
        tval = hybrid_model.time_symbolic_models[1].tsteps[aug[2]]
        @test 0.0 <= tval <= 1.0
        xval = SY.TimedHybridAutomata.get_concrete_state(hybrid_model, idx)[1][1]
        @test 0.0 <= xval <= 0.5
    end
end

@testset "TimedHybridAutomata - simple, Time not taken into account" begin
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

    growth_bounds = SVector(SMatrix{1, 1}(0.5), SMatrix{1, 1}(0.8))

    param_discretisation = [
        (0.5, 0.5, 0.5),  # mode 1
        (0.5, 0.3, 0.5),   # mode 2
    ]

    hybrid_model = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        growth_bounds,
        param_discretisation,
    )

    # ===================== TESTS =====================

    for aug_state in hybrid_model.int2aug_state
        # aug_state = (state_symbol, time_symbol, mode_id)
        @test aug_state[2] == 1
    end

    for tm in hybrid_model.time_symbolic_models
        @test length(tm.tsteps) == 1
        @test tm.is_active == false
    end

    for tm in hybrid_model.time_symbolic_models
        tmin, tmax = tm.tsteps[1], tm.tsteps[1]
        for tval in [tmin, tmax, (tmin + tmax)/2]
            @test SY.TimedHybridAutomata.find_time_index(tm, tval) == 1
        end
    end

    for tm in hybrid_model.time_symbolic_models
        tmin, tmax = tm.tsteps[1], tm.tsteps[1]
        @test SY.TimedHybridAutomata.get_time_indices_from_interval(tm, [tmin, tmax]) == [1]
    end

    @test SY.TimedHybridAutomata.get_n_state(hybrid_model) ==
          length(hybrid_model.int2aug_state)
    @test collect(SY.TimedHybridAutomata.enum_states(hybrid_model)) ==
          collect(1:SY.TimedHybridAutomata.get_n_state(hybrid_model))

    @test SY.TimedHybridAutomata.get_n_input(hybrid_model) ==
          hybrid_model.global_input_map.total_inputs

    for k in 1:length(hybrid_model.symmodels)
        expected = collect(
            1:SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(
                hybrid_model.symmodels[k],
            ),
        )
        actual = collect(SY.TimedHybridAutomata.enum_inputs(hybrid_model, k))
        @test actual == expected
    end

    for state in 1:SY.TimedHybridAutomata.get_n_state(hybrid_model)
        aug = hybrid_model.int2aug_state[state]
        concrete = SY.TimedHybridAutomata.get_concrete_state(hybrid_model, state)
        idx = SY.TimedHybridAutomata.get_abstract_state(hybrid_model, concrete)
        @test idx == state
    end

    gim = hybrid_model.global_input_map
    for k in 1:length(hybrid_model.symmodels)
        n_inputs =
            SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(hybrid_model.symmodels[k])
        for local_input_id in 1:n_inputs
            global_input_id =
                SY.TimedHybridAutomata.get_global_input_id(gim, k, local_input_id)

            u_concrete =
                SY.TimedHybridAutomata.get_concrete_input(hybrid_model, global_input_id, k)
            global_input_id2 =
                SY.TimedHybridAutomata.get_abstract_input(hybrid_model, u_concrete, k)
            @test global_input_id == global_input_id2
        end
    end

    for tr in 1:gim.switching_inputs
        gid = SY.TimedHybridAutomata.get_switching_global_id(gim, tr)
        @test isnothing(SY.TimedHybridAutomata.get_concrete_input(hybrid_model, gid, 1))
    end

    if hasproperty(gim, :switch_labels)
        @test length(gim.switch_labels) == gim.switching_inputs
        for label in gim.switch_labels
            @test occursin("SWITCH", label)
            @test occursin("->", label)
        end
    end
end

@testset "TimedHybridAutomata - 3 modes, transitions complètes, reset maps variées" begin
    #
    # Test 3: 3-mode hybrid system (2D state + 1D time)
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

    growth_bounds = SVector(
        SMatrix{2, 2}(abs.(A_A)),
        SMatrix{2, 2}(abs.(A_B)),
        SMatrix{2, 2}(abs.(A_C)),
    )
    param_discretisation = [
        (0.5, 0.5, 0.5),  # mode A
        (0.5, 0.5, 0.5),  # mode B
        (0.5, 0.5, 0.5),   # mode C
    ]

    hybrid_model = SY.TimedHybridAutomata.Build_Timed_Hybrid_Automaton(
        hs,
        growth_bounds,
        param_discretisation,
    )
    # ===================== TESTS =====================

    @test length(hybrid_model.symmodels) == 3
    @test length(hybrid_model.time_symbolic_models) == 3
    @test length(hybrid_model.int2aug_state) == length(hybrid_model.aug_state2int)

    gim = hybrid_model.global_input_map
    @test gim.switching_inputs == 6
    @test length(gim.continuous_to_global) == 25 + 25 + 7*7
    @test length(gim.switching_to_global) == 6

    for mode in 1:3
        for local_input in 1:3
            gid = SY.TimedHybridAutomata.get_global_input_id(gim, mode, local_input)
            @test gid > 0
            typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, gid)
            @test typ == :continuous
            @test info[1] == mode
        end
    end
    for tr in 1:6
        gid = SY.TimedHybridAutomata.get_switching_global_id(gim, tr)
        @test gid > 0
        typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, gid)
        @test typ == :switching
        @test info == tr
    end

    mode_n_input = [25, 7*7, 25]
    for (idxsm, sm) in enumerate(hybrid_model.symmodels)
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_state(sm) > 0
        @test SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(sm) ==
              mode_n_input[idxsm]
    end
    for tm in hybrid_model.time_symbolic_models
        @test length(tm.tsteps) > 0
        @test tm.tsteps[1] <= tm.tsteps[end]
    end

    for i in 1:length(hybrid_model.int2aug_state)
        @test haskey(hybrid_model.aug_state2int, hybrid_model.int2aug_state[i])
    end

    for mode in 1:3
        sm = hybrid_model.symmodels[mode]
        tm = hybrid_model.time_symbolic_models[mode]
        x0 = [0.0, 0.0]
        idx = SY.TimedHybridAutomata.find_symbolic_state(sm, x0)
        @test idx > 0
        t0 = tm.tsteps[1]
        tidx = SY.TimedHybridAutomata.find_time_index(tm, t0)
        @test tidx == 1
        idx_bad = SY.TimedHybridAutomata.find_symbolic_state(sm, [100.0, 100.0])
        @test idx_bad == 0
        guard = guards[mode]
        sp = SY.TimedHybridAutomata.extract_spatial_part(guard)
        tp = SY.TimedHybridAutomata.extract_temporal_part(guard)
        @test sp.lb == guard.lb[1:(end - 1)]
        @test tp == [guard.lb[end], guard.ub[end]]
        indices = SY.TimedHybridAutomata.get_time_indices_from_interval(tm, tp)
        @test all(i >= 1 && i <= length(tm.tsteps) for i in indices)
    end

    for mode in 1:3
        sm = hybrid_model.symmodels[mode]
        @test SY.TimedHybridAutomata.find_symbolic_state(sm, [999.0, -999.0]) == 0
    end

    for gid in 1:gim.total_inputs
        typ, info = SY.TimedHybridAutomata.get_local_input_info(gim, gid)
        if typ == :continuous
            @test 1 <= info[1] <= 27
            @test 1 <= info[2] <= 99
        elseif typ == :switching
            @test 1 <= info <= 6
        else
            @test typ == :invalid
        end
    end

    for i in 1:length(hybrid_model.int2aug_state)
        aug = hybrid_model.int2aug_state[i]
        mode = aug[3]
        t_idx = aug[2]
        @test 1 <= mode <= 3
        @test 1 <= t_idx <= length(hybrid_model.time_symbolic_models[mode].tsteps)
    end

    # --- Tests des accessors et des labels de commutation ---

    @test SY.TimedHybridAutomata.get_n_state(hybrid_model) ==
          length(hybrid_model.int2aug_state)
    @test SY.TimedHybridAutomata.get_n_input(hybrid_model) ==
          hybrid_model.global_input_map.total_inputs

    @test collect(SY.TimedHybridAutomata.enum_states(hybrid_model)) ==
          collect(1:SY.TimedHybridAutomata.get_n_state(hybrid_model))

    for k in 1:length(hybrid_model.symmodels)
        expected = collect(
            1:SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(
                hybrid_model.symmodels[k],
            ),
        )
        actual = collect(SY.TimedHybridAutomata.enum_inputs(hybrid_model, k))
        @test actual == expected
    end

    for state in 1:SY.TimedHybridAutomata.get_n_state(hybrid_model)
        aug = hybrid_model.int2aug_state[state]
        concrete = SY.TimedHybridAutomata.get_concrete_state(hybrid_model, state)
        idx = SY.TimedHybridAutomata.get_abstract_state(hybrid_model, concrete)
        @test idx == state
    end

    gim = hybrid_model.global_input_map
    for k in 1:length(hybrid_model.symmodels)
        n_inputs =
            SY.TimedHybridAutomata.Dionysos.Symbolic.get_n_input(hybrid_model.symmodels[k])
        for local_input_id in 1:n_inputs
            global_input_id =
                SY.TimedHybridAutomata.get_global_input_id(gim, k, local_input_id)
            u_concrete =
                SY.TimedHybridAutomata.get_concrete_input(hybrid_model, global_input_id, k)
            global_input_id2 =
                SY.TimedHybridAutomata.get_abstract_input(hybrid_model, u_concrete, k)
            @test global_input_id == global_input_id2
        end
    end

    for tr in 1:gim.switching_inputs
        gid = SY.TimedHybridAutomata.get_switching_global_id(gim, tr)
        @test isnothing(SY.TimedHybridAutomata.get_concrete_input(hybrid_model, gid, 1))
    end

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
end
end
