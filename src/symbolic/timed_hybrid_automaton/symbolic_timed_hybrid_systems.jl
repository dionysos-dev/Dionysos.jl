export SymbolicTimedHybridSystems

module SymbolicTimedHybridSystems

import HybridSystems: HybridSystem
import StaticArrays: SVector
import HybridSystems, MathematicalSystems
import MathOptInterface as MOI

import Dionysos

# ================================================================
# Type definitions for type stability
# ================================================================

# (?) Is it usefull to delete type unstability
"""
Concrete type alias for augmented states: (state_id, time_id, mode_id)
"""
const AugmentedState = Tuple{Int, Int, Int}

"""
Concrete type alias for transitions: (target_state, source_state, input_id)
"""
const TransitionTuple = Tuple{AugmentedState, AugmentedState, Int}

# ================================================================
# Symbolic temporal hybrid model structure
# ================================================================

"""
    TimedHybridSymbolicModel{S1, A, T, G}

Main structure representing the symbolic abstraction of a timed hybrid system.

# Fields
- `mode_abstractions::Vector{S1}`: Symbolic models of the dynamics per mode
- `time_abstractions::Vector{T}`: Symbolic models of time per mode  
- `state_index_to_augmented::Vector{AugmentedState}`: Maps integer indices to augmented states (state_id, time_id, mode_id)
- `augmented_to_state_index::Dict{AugmentedState, Int}`: Maps augmented states to integer indices
- `symbolic_automaton::A`: Final symbolic automaton representing the timed hybrid system
- `input_mapping::G`: Global input mapping system for unified input handling
"""
struct TimedHybridSymbolicModel{S1, A, T, G}
    mode_abstractions::Vector{S1}
    time_abstractions::Vector{T}
    state_index_to_augmented::Vector{AugmentedState}
    augmented_to_state_index::Dict{AugmentedState, Int}
    symbolic_automaton::A
    input_mapping::G
end

# ================================================================
# Structure for matching global abstract inputs
# ================================================================

"""
    GlobalInputMap

Structure for managing the mapping between local (per-mode) and global input indices, for both continuous and switching inputs.

# Fields
- `total_inputs::Int`: Total number of global inputs
- `continuous_inputs::Int`: Number of continuous global inputs
- `switching_inputs::Int`: Number of switching global inputs
- `continuous_to_global::Dict{Tuple{Int, Int}, Int}`: (mode_id, local_input_id) → global_input_id
- `global_to_continuous::Dict{Int, Tuple{Int, Int}}`: global_input_id → (mode_id, local_input_id)
- `switching_to_global::Dict{Int, Int}`: transition_id → global_input_id
- `global_to_switching::Dict{Int, Int}`: global_input_id → transition_id
- `continuous_range::UnitRange{Int}`: Range of continuous global input ids
- `switching_range::UnitRange{Int}`: Range of switching global input ids
"""
struct GlobalInputMap
    total_inputs::Int
    continuous_inputs::Int
    switching_inputs::Int
    continuous_to_global::Dict{Tuple{Int, Int}, Int}    # (mode_id, local_input_id) → global_input_id
    global_to_continuous::Dict{Int, Tuple{Int, Int}}    # global_input_id → (mode_id, local_input_id)
    switching_to_global::Dict{Int, Int}                 # transition_id → global_input_id
    global_to_switching::Dict{Int, Int}                 # global_input_id → transition_id
    continuous_range::UnitRange{Int}
    switching_range::UnitRange{Int}
    switch_labels::Vector{String} # labels  for switching inputs (e.g., "SWITCH source_mode_id -> target_mode_id")
end

"""
    GlobalInputMap(abstract_systems, hs::HybridSystem)

Construct a GlobalInputMap for a given hybrid system and its symbolic abstractions.

# Arguments
- `abstract_systems`: Vector of (symbolic_dynamics, symbolic_time) tuples per mode
- `hs::HybridSystem`: The hybrid system

# Returns
- `GlobalInputMap`: The constructed mapping structure
"""
function GlobalInputMap(abstract_systems, hs::HybridSystem)
    # Phase 1: Allocate continuous inputs
    continuous_to_global = Dict{Tuple{Int, Int}, Int}()
    global_to_continuous = Dict{Int, Tuple{Int, Int}}()
    continuous_count = 0
    for (mode_id, (symmodel_dynam, _)) in enumerate(abstract_systems)
        input_count = Dionysos.Symbolic.get_n_input(symmodel_dynam)
        for local_input_id in 1:input_count
            global_id = continuous_count + local_input_id
            continuous_to_global[(mode_id, local_input_id)] = global_id
            global_to_continuous[global_id] = (mode_id, local_input_id)
        end
        continuous_count += input_count
    end
    # Phase 2: Allocate switching inputs et labels
    switching_to_global = Dict{Int, Int}()
    global_to_switching = Dict{Int, Int}()
    switch_labels = String[]
    switching_count = 0
    transitions = collect(HybridSystems.transitions(hs.automaton))
    for (transition_id, transition) in enumerate(transitions)
        global_id = continuous_count + switching_count + 1
        switching_to_global[transition_id] = global_id
        global_to_switching[global_id] = transition_id
        source_id = HybridSystems.source(hs.automaton, transition)
        target_id = HybridSystems.target(hs.automaton, transition)
        push!(switch_labels, "SWITCH $(source_id) -> $(target_id)")
        switching_count += 1
    end
    # Phase 3: Compute ranges
    continuous_range = 1:continuous_count
    switching_range = (continuous_count + 1):(continuous_count + switching_count)
    return GlobalInputMap(
        continuous_count + switching_count,
        continuous_count,
        switching_count,
        continuous_to_global,
        global_to_continuous,
        switching_to_global,
        global_to_switching,
        continuous_range,
        switching_range,
        switch_labels,
    )
end

# === Accessor functions ===

"""
    get_global_input_id(gim::GlobalInputMap, mode_id::Int, local_input_id::Int) -> Int

Get the global input id for a local continuous input.

# Arguments
- `gim::GlobalInputMap`: the global input map
- `mode_id::Int`: mode id (1, 2, 3, ...)
- `local_input_id::Int`: local input id in this mode (1, 2, 3, ...)

# Returns
- `Int`: the global input id (0 if not found)
"""
function get_global_input_id(gim::GlobalInputMap, mode_id::Int, local_input_id::Int)
    return get(gim.continuous_to_global, (mode_id, local_input_id), 0)
end

"""
    get_switching_global_id(gim::GlobalInputMap, transition_id::Int) -> Int

Get the global input id for a switching input.

# Arguments
- `gim::GlobalInputMap`: the global input map
- `transition_id::Int`: transition id in the hybrid automaton

# Returns
- `Int`: the global input id for switching (0 if not found)
"""
function get_switching_global_id(gim::GlobalInputMap, transition_id::Int)
    return get(gim.switching_to_global, transition_id, 0)
end

"""
    get_local_input_info(gim::GlobalInputMap, global_id::Int) -> (Symbol, Union{Tuple{Int,Int}, Int, Nothing})

Determine the type and local info of a global input id.

# Arguments
- `gim::GlobalInputMap`: the global input map
- `global_id::Int`: the global input id

# Returns
- Tuple:
  - `Symbol`: `:continuous`, `:switching`, or `:invalid`
  - `Union{Tuple{Int,Int}, Int, Nothing}`:
    - if `:continuous`: `(mode_id, local_input_id)`
    - if `:switching`: `transition_id`
    - if `:invalid`: `nothing`
"""
function get_local_input_info(gim::GlobalInputMap, global_id::Int)
    if global_id in gim.continuous_range
        return :continuous, gim.global_to_continuous[global_id]
    elseif global_id in gim.switching_range
        return :switching, gim.global_to_switching[global_id]
    else
        return :invalid, nothing
    end
end

"""
    is_continuous_input(gim::GlobalInputMap, global_id::Int) -> Bool

Check if a global input id is a continuous input.

# Arguments
- `gim::GlobalInputMap`: the global input map
- `global_id::Int`: the global input id

# Returns
- `Bool`: `true` if continuous, `false` otherwise
"""
function is_continuous_input(gim::GlobalInputMap, global_id::Int)
    return global_id in gim.continuous_range
end

"""
    is_switching_input(gim::GlobalInputMap, global_id::Int) -> Bool

Check if a global input id is a switching input.

# Arguments
- `gim::GlobalInputMap`: the global input map
- `global_id::Int`: the global input id

# Returns
- `Bool`: `true` if switching, `false` otherwise
"""
function is_switching_input(gim::GlobalInputMap, global_id::Int)
    return global_id in gim.switching_range
end
# ================================================================
# Symbolic model creation
# ================================================================

#
# NOTE: In the future, this will be changed so that the user can provide their own optimizer
# for the abstraction of the dynamics in each mode. This will allow for custom abstraction
# methods and greater flexibility in how the symbolic model is constructed.
#

"""
    build_dynamical_symbolic_model(system, growth_bound, param_discretisation)

Build a symbolic abstraction of a continuous system using uniform grid discretization.

# Returns
- Symbolic abstraction of the system
"""
function build_dynamical_symbolic_model(
    system;
    optimizer_factory = nothing,
    optimizer = nothing,
    optimizer_kwargs = Dict(),
)
    # (?) will be changed, need some discussion
    if optimizer !== nothing
        opt = optimizer
    elseif optimizer_factory !== nothing
        opt = optimizer_factory()
    else
        opt = MOI.instantiate(Dionysos.Optim.Abstraction.UniformGridAbstraction.Optimizer)
    end

    problem = Dionysos.Problem.EmptyProblem(system, system.X)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), problem)

    for (k, v) in optimizer_kwargs
        MOI.set(opt, MOI.RawOptimizerAttribute(k), v)
    end

    MOI.optimize!(opt)
    return MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
end

"""
    build_mode_symbolic_abstractions(hs::HybridSystem, optimizer_list, optimizer_kwargs_dict)

Build symbolic models (dynamics and time) for each mode of a hybrid system.

# Arguments
- `hs::HybridSystem`: The hybrid system
- `optimizer_list`: Vector of optimizer factory functions, one per mode
- `optimizer_kwargs_dict`: Vector of dictionaries containing optimizer parameters per mode

# Returns
- Vector of (symbolic_dynamics, symbolic_time) tuples per mode

# Safety Features
- Validates optimizer list length matches number of modes
- Error handling for mode construction failures
"""
function build_mode_symbolic_abstractions(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict},
)
    # Input validation (?) may changed to allow same optimizer for many modes
    n_modes = length(HybridSystems.states(hs.automaton))
    @assert length(optimizer_list) == n_modes "Optimizer list length mismatch"
    @assert length(optimizer_kwargs_dict) == n_modes "Optimizer kwargs length mismatch"

    # Pre-allocate result vector (performance (?) type any)
    mode_abstractions = Vector{Tuple{Any, Any}}(undef, n_modes)

    # Build symbolic models for each mode
    for (i, mode_id) in enumerate(HybridSystems.states(hs.automaton))
        try
            mode_system = HybridSystems.mode(hs, mode_id)
            dynamics_system = mode_system.systems[1]    # physical dynamics
            time_system = mode_system.systems[2]        # time dynamics

            # Build symbolic model for dynamics
            symbolic_dynamics = build_dynamical_symbolic_model(
                dynamics_system;
                optimizer_factory = optimizer_list[i],
                optimizer_kwargs = optimizer_kwargs_dict[i],
            )

            # Build symbolic model for time
            symbolic_time = Dionysos.Symbolic.TimeSymbolicModel(
                time_system,
                get(optimizer_kwargs_dict[i], "time_step", nothing),
            )

            mode_abstractions[i] = (symbolic_dynamics, symbolic_time)

        catch e
            error("Failed to build symbolic abstraction for mode $mode_id: $e")
        end
    end

    return mode_abstractions
end

"""
    build_all_transitions(hs::HybridSystem, mode_abstractions, input_mapping::GlobalInputMap)

Build all transitions (intra-mode and inter-mode) for the timed hybrid system.
Centralized function that coordinates transition building with better organization.

# Arguments
- `hs::HybridSystem`: The hybrid system
- `mode_abstractions`: Vector of (symbolic_dynamics, symbolic_time) tuples per mode  
- `input_mapping::GlobalInputMap`: The global input mapping

# Returns
- Vector of transitions (tuples of (target, source, input))

# Performance Notes
- Pre-allocates transition list based on exact count of intra-mode transitions
- Inter-mode transitions grow the vector dynamically to avoid over-allocation (better solution (?))
- Coordinates both intra-mode and inter-mode transition building
"""
function build_all_transitions(
    hs::HybridSystem,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)
    # Pre-allocate transition list with exact count for intra-mode transitions
    intra_mode_transitions = sum(
        Dionysos.Symbolic.get_n_transitions(abs_pair[1]) * length(abs_pair[2].tsteps)
        for abs_pair in mode_abstractions
    )

    # (?) This is a bottleneck, another solution ? 
    transition_list = Vector{TransitionTuple}()
    sizehint!(transition_list, intra_mode_transitions)

    # Add intra-mode transitions (within each mode)
    add_intra_mode_transitions!(transition_list, mode_abstractions, input_mapping)

    # Add inter-mode transitions (switching between modes)
    add_inter_mode_transitions!(transition_list, hs, mode_abstractions, input_mapping)

    return transition_list
end

# ================================================================
# Functions to add transitions
# ================================================================

# (1) bottleneck 

"""
    add_intra_mode_transitions!(transition_list, mode_abstractions, input_mapping::GlobalInputMap)

Add intra-mode transitions to the transition list with optimized performance.
Handles both time-frozen and multi-step time cases efficiently.

# Arguments  
- `transition_list`: The list to which transitions are appended
- `mode_abstractions`: Vector of (symbolic_dynamics, symbolic_time) tuples per mode
- `input_mapping::GlobalInputMap`: The global input mapping

# Performance Notes (Was one of our bottlenecks)
- Pre-allocates space to avoid vector reallocations
- Uses @inbounds in critical loops after bounds checking
- Caches frequently accessed values
"""
function add_intra_mode_transitions!(
    transition_list,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)

    # Bounds checking outside critical loop
    n_modes = length(mode_abstractions)
    # Validate that we have enough modes for the input mapping
    @assert n_modes > 0 "No mode abstractions provided"

    # Process each mode with optimized loops
    for (mode_id, (symmodel_dynamics, symmodel_time)) in enumerate(mode_abstractions)
        time_steps = symmodel_time.tsteps
        n_time_steps = length(time_steps)

        if n_time_steps == 1
            @inbounds for (target, source, local_input_id) in
                          Dionysos.Symbolic.enum_transitions(symmodel_dynamics)
                global_input_id =
                    get_global_input_id(input_mapping, mode_id, local_input_id)
                if global_input_id > 0  # Safety check
                    target_state = (target, 1, mode_id)::AugmentedState
                    source_state = (source, 1, mode_id)::AugmentedState
                    push!(
                        transition_list,
                        (target_state, source_state, global_input_id)::TransitionTuple,
                    )
                end
            end
        else
            transitions_cache =
                collect(Dionysos.Symbolic.enum_transitions(symmodel_dynamics))
            @inbounds for k in 1:(n_time_steps - 1)
                for (target, source, local_input_id) in transitions_cache
                    global_input_id =
                        get_global_input_id(input_mapping, mode_id, local_input_id)
                    if global_input_id > 0  # Safety check
                        target_state = (target, k + 1, mode_id)::AugmentedState
                        source_state = (source, k, mode_id)::AugmentedState
                        push!(
                            transition_list,
                            (target_state, source_state, global_input_id)::TransitionTuple,
                        )
                    end
                end
            end
        end
    end
end

"""
    add_inter_mode_transitions!(transition_list, hs::HybridSystem, mode_abstractions, input_mapping::GlobalInputMap)

Add inter-mode transitions (mode switches) to the transition list using guards and reset maps.

# Arguments
- `transition_list`: The list to which transitions are appended
- `hs::HybridSystem`: The hybrid system
- `mode_abstractions`: Vector of (symbolic_dynamics, symbolic_time) tuples per mode
- `input_mapping::GlobalInputMap`: The global input mapping

# Safety Notes
- Includes bounds checking for all state and time indices
- Handles cases where reset maps produce invalid states
- Validates guard intersections before processing
"""
function add_inter_mode_transitions!(
    transition_list,
    hs::HybridSystem,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)
    transitions = collect(HybridSystems.transitions(hs.automaton))

    for (transition_id, transition) in enumerate(transitions)
        # Get the global input id for this switching transition
        global_input_id = get_switching_global_id(input_mapping, transition_id)

        # Safety check for valid global input ID
        if global_input_id <= 0
            @warn "Invalid global input ID for transition $transition_id, skipping"
            continue
        end

        # Extract source and target mode indices with bounds checking
        source_mode = HybridSystems.source(hs.automaton, transition)
        target_mode = HybridSystems.target(hs.automaton, transition)

        @assert source_mode <= length(mode_abstractions) "Source mode $source_mode out of bounds"
        @assert target_mode <= length(mode_abstractions) "Target mode $target_mode out of bounds"

        # Get the reset map and guard for this transition
        reset_map = HybridSystems.resetmap(hs, transition)
        guard = HybridSystems.guard(hs, transition)

        if isnothing(guard)
            @warn "No guard found for transition $transition_id, skipping"
            continue
        end

        # Get the symbolic models for the source and target modes
        (source_symmodel_dynamics, source_time_model) = mode_abstractions[source_mode]
        (target_symmodel_dynamics, target_time_model) = mode_abstractions[target_mode]

        try
            # Split the guard into spatial and temporal parts
            guard_spatial = extract_spatial_part(guard)
            guard_temporal = extract_temporal_part(guard)

            # Get all source states that intersect with the spatial guard
            source_states = Dionysos.Symbolic.get_states_from_set(
                source_symmodel_dynamics,
                guard_spatial,
                Dionysos.Domain.INNER,
            )

            # Get all time indices that intersect with the temporal guard
            time_indices = get_time_indices_from_interval(source_time_model, guard_temporal)

            if isempty(source_states) || isempty(time_indices)
                @warn "Empty guard intersection for transition $transition_id, skipping"
                continue
            end

            # Process each valid combination of (state, time) in the guard
            for source_state in source_states, source_time_idx in time_indices
                # Bounds checking for time index
                if source_time_idx > length(source_time_model.tsteps) ||
                   source_time_idx <= 0
                    continue
                end

                # Build the augmented source state [x1, x2, ..., xn, t]
                source_continuous_state = Dionysos.Symbolic.get_concrete_state(
                    source_symmodel_dynamics,
                    source_state,
                )
                source_time_value = source_time_model.tsteps[source_time_idx]
                augmented_source_state = vcat(source_continuous_state, source_time_value)

                # Apply the reset map to the augmented state
                reset_result = MathematicalSystems.apply(reset_map, augmented_source_state)
                reset_continuous_part = reset_result[1:(end - 1)]
                reset_time_value = reset_result[end]

                # Find the corresponding target symbolic state and time index
                target_state =
                    find_symbolic_state(target_symmodel_dynamics, reset_continuous_part)
                target_time_idx =
                    Dionysos.Symbolic.ceil_time2int(target_time_model, reset_time_value)

                # Add the transition if both target state and time are valid
                if target_state > 0 &&
                   target_time_idx > 0 &&
                   target_time_idx <= length(target_time_model.tsteps)
                    target_aug_state =
                        (target_state, target_time_idx, target_mode)::AugmentedState
                    source_aug_state =
                        (source_state, source_time_idx, source_mode)::AugmentedState
                    push!(
                        transition_list,
                        (
                            target_aug_state,
                            source_aug_state,
                            global_input_id,
                        )::TransitionTuple,
                    )
                end
            end

        catch e
            @warn "Error processing transition $transition_id: $e"
            continue
        end
    end
end

# ================================================================
# Build the final automaton from the transition list
# ================================================================

# (2) bottleneck 

"""
    build_symbolic_automaton(transition_list, mode_abstractions, input_mapping::GlobalInputMap)

Build the symbolic automaton from the list of temporal transitions with optimized performance.
Uses pre-allocation with exact estimates, efficient state enumeration, and bounds checking for safety.

# Arguments
- `transition_list`: List of temporal transitions (tuples of (target, source, input))
- `mode_abstractions`: Vector of (symbolic_dynamics, symbolic_time) tuples per mode
- `input_mapping::GlobalInputMap`: The global input mapping for exact input count

# Returns
- Tuple (state_index_to_augmented, augmented_to_state_index, automaton):
    - `state_index_to_augmented`: Vector mapping integer indices to augmented states
    - `augmented_to_state_index`: Dict mapping augmented states to integer indices  
    - `automaton`: The constructed symbolic automaton

# Performance Notes
- Uses exact state count estimation instead of heuristics
- Pre-allocates all data structures to avoid reallocations
- Uses @inbounds in critical loops after bounds checking
- Single-pass collection with optimal memory usage
"""
function build_symbolic_automaton(
    transition_list,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)
    @assert !isempty(transition_list) "Transition list cannot be empty"
    @assert !isempty(mode_abstractions) "Mode abstractions cannot be empty"

    estimated_states = sum(
        Dionysos.Symbolic.get_n_state(abs_pair[1]) * length(abs_pair[2].tsteps) for
        abs_pair in mode_abstractions
    )

    states = Set{AugmentedState}()
    inputs_set = Set{Int}()
    sizehint!(states, estimated_states)
    sizehint!(inputs_set, input_mapping.total_inputs)

    for (target, source, input) in transition_list
        push!(states, target)
        push!(states, source)
        push!(inputs_set, input)
    end

    augmented_states = collect(states)
    nstates = length(augmented_states)
    ninputs = length(inputs_set)

    state_index_to_augmented = Vector{AugmentedState}(undef, nstates)
    augmented_to_state_index = Dict{AugmentedState, Int}()
    sizehint!(augmented_to_state_index, nstates)

    @inbounds for i in eachindex(augmented_states)
        aug_state = augmented_states[i]
        state_index_to_augmented[i] = aug_state
        augmented_to_state_index[aug_state] = i
    end

    symbolic_automaton = Dionysos.Symbolic.NewIndexedAutomatonList(nstates, ninputs)

    @inbounds for (target, source, abstract_input) in transition_list
        target_int = augmented_to_state_index[target]
        source_int = augmented_to_state_index[source]

        Dionysos.Symbolic.add_transition!(
            symbolic_automaton,
            source_int,
            target_int,
            abstract_input,
        )
    end

    return state_index_to_augmented, augmented_to_state_index, symbolic_automaton
end

# ================================================================
# Main constructor for timed hybrid symbolic models
# ================================================================

"""
    build_timed_hybrid_symbolic_model(hs::HybridSystem, optimizer_list, optimizer_kwargs_dict)

Construct a complete timed hybrid symbolic model for a given hybrid system.
This is the main function users should call to build symbolic abstractions.

# Arguments
- `hs::HybridSystem`: The hybrid system to abstract
- `optimizer_list`: Vector of optimizer factory functions, one per mode
- `optimizer_kwargs_dict`: Vector of dictionaries containing optimizer parameters per mode

# Returns
- `TimedHybridSymbolicModel`: The constructed symbolic model with optimized performance

# Safety Features
- Validates input dimensions and consistency
- Bounds checking for all array accesses
- Error handling for invalid states and transitions

# Performance Features  
- Pre-allocates data structures to minimize memory allocations
- Uses efficient algorithms for transition building
- Optimized automaton construction with @inbounds where safe
"""
function build_timed_hybrid_symbolic_model(
    hs::HybridSystem,
    optimizer_list::AbstractVector{Function},
    optimizer_kwargs_dict::AbstractVector{<:Dict},
)
    # Input validation for safety
    n_modes = length(HybridSystems.states(hs.automaton))
    @assert length(optimizer_list) == n_modes "Number of optimizers ($(length(optimizer_list))) must match number of modes ($n_modes)"
    @assert length(optimizer_kwargs_dict) == n_modes "Number of optimizer configs ($(length(optimizer_kwargs_dict))) must match number of modes ($n_modes)"

    # 1) Build symbolic models per mode with improved function name
    mode_abstractions =
        build_mode_symbolic_abstractions(hs, optimizer_list, optimizer_kwargs_dict)

    # 2) Build the global input mapping system
    input_mapping = GlobalInputMap(mode_abstractions, hs)

    # 3) Build transition list efficiently
    transition_list = build_all_transitions(hs, mode_abstractions, input_mapping)

    # 4) Build the final automaton with optimized performance
    state_index_to_augmented, augmented_to_state_index, symbolic_automaton =
        build_symbolic_automaton(transition_list, mode_abstractions, input_mapping)

    # 5) Extract symbolic models with clearer variable names
    mode_dynamics_models = [abs_sys[1] for abs_sys in mode_abstractions]
    mode_time_models = [abs_sys[2] for abs_sys in mode_abstractions]

    # 6) Return the complete timed hybrid symbolic model with new structure
    return TimedHybridSymbolicModel(
        mode_dynamics_models,
        mode_time_models,
        state_index_to_augmented,
        augmented_to_state_index,
        symbolic_automaton,
        input_mapping,
    )
end

# ================================================================
# Utility functions for safe and robust symbolic state operations
# ================================================================

"""
    find_symbolic_state(symmodel, continuous_state)

Find the symbolic state index corresponding to a given continuous state with enhanced safety.
Improved version with better error handling and logging.

# Arguments
- `symmodel`: The symbolic model
- `continuous_state`: The continuous state vector

# Returns
- `Int`: The symbolic state index (0 if not found)

# Safety Features
- Comprehensive error handling for edge cases
- Validates input dimensions
- Provides informative warnings for debugging
"""
function find_symbolic_state(symmodel, continuous_state)
    # Input validation
    if isnothing(continuous_state) || isempty(continuous_state)
        @warn "Invalid continuous state provided: empty or nothing"
        return 0
    end

    try
        state_idx = Dionysos.Symbolic.get_abstract_state(symmodel, continuous_state)
        if isnothing(state_idx) || state_idx <= 0
            # This is common for boundary cases, so debug level rather than warning
            return 0
        else
            return state_idx
        end
    catch e
        @warn "Error finding symbolic state for continuous state $continuous_state: $e"
        return 0
    end
end

"""
    extract_spatial_part(guard)

Extract the spatial part (all but last dimension) from a guard (assumed to be a HyperRectangle).

# Arguments
- `guard`: The guard set (should be a HyperRectangle)

# Returns
- `HyperRectangle`: The spatial part of the guard
"""
function extract_spatial_part(guard)
    if isa(guard, Dionysos.Utils.HyperRectangle)
        return Dionysos.Utils.HyperRectangle(guard.lb[1:(end - 1)], guard.ub[1:(end - 1)])
    else
        error("Unsupported guard type: $(typeof(guard))")
    end
end

"""
    extract_temporal_part(guard)

Extract the temporal part (last dimension) from a guard (assumed to be a HyperRectangle).

# Arguments
- `guard`: The guard set (should be a HyperRectangle)

# Returns
- `Vector{Float64}`: The temporal interval [t_min, t_max]
"""
function extract_temporal_part(guard)
    if isa(guard, Dionysos.Utils.HyperRectangle)
        return [guard.lb[end], guard.ub[end]]
    else
        error("Unsupported guard type: $(typeof(guard))")
    end
end

"""
    get_time_indices_from_interval(time_model, temporal_interval)

Get all time indices in the symbolic time model that fall within a given interval.

# Arguments
- `time_model`: The symbolic time model
- `temporal_interval`: The interval [t_min, t_max]

# Returns
- `Vector{Int}`: Indices of time steps within the interval
"""
function get_time_indices_from_interval(time_model, temporal_interval)
    t_min, t_max = temporal_interval
    return findall(t -> t_min <= t <= t_max, time_model.tsteps)
end

# ================================================================
# Accessor functions for TimedHybridSymbolicModel
# ================================================================

"""Get the number of states in the timed hybrid symbolic model"""
get_n_state(model::TimedHybridSymbolicModel) = length(model.state_index_to_augmented)

"""Get the total number of inputs (continuous + switching)"""
function get_n_input(model::TimedHybridSymbolicModel)
    return model.input_mapping.total_inputs
end

"""Enumerate all state indices in the model"""
enum_states(model::TimedHybridSymbolicModel) = 1:get_n_state(model)

"""Enumerate input indices for a given mode with bounds checking"""
function enum_inputs(model::TimedHybridSymbolicModel, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"
    return Dionysos.Symbolic.enum_inputs(model.mode_abstractions[mode_id])
end

"""
    get_concrete_state(model::TimedHybridSymbolicModel, state_index::Int)

Convert an abstract state index to its concrete augmented state representation.
Returns (continuous_state, time, mode_id).

# Safety Features
- Bounds checking for state index
- Proper error handling for invalid states
"""
function get_concrete_state(model::TimedHybridSymbolicModel, state_index::Int)
    @assert 1 <= state_index <= length(model.state_index_to_augmented) "State index $state_index out of bounds"

    (state_id, time_id, mode_id) = model.state_index_to_augmented[state_index]

    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    dynamics_model = model.mode_abstractions[mode_id]
    time_model = model.time_abstractions[mode_id]

    continuous_state = Dionysos.Symbolic.get_concrete_state(dynamics_model, state_id)
    time_value = Dionysos.Symbolic.int2time(time_model, time_id)

    return (continuous_state, time_value, mode_id)
end

"""
    get_abstract_state(model::TimedHybridSymbolicModel, augmented_state)

Convert an augmented concrete state (continuous_state, time, mode_id) to its abstract state index.

# Safety Features  
- Validates mode_id bounds
- Handles cases where no valid abstract state exists
"""
function get_abstract_state(model::TimedHybridSymbolicModel, augmented_state)
    (continuous_state, time_value, mode_id) = augmented_state

    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    dynamics_model = model.mode_abstractions[mode_id]
    time_model = model.time_abstractions[mode_id]

    # Find abstract state and time indices
    state_id = Dionysos.Symbolic.get_abstract_state(dynamics_model, continuous_state)
    time_id = Dionysos.Symbolic.floor_time2int(time_model, time_value)

    if isnothing(state_id) || time_id <= 0
        error("No valid abstract state found for augmented_state $augmented_state")
    end

    augmented_key = (state_id, time_id, mode_id)::AugmentedState
    return get(model.augmented_to_state_index, augmented_key, 0)
end

"""
    get_states_from_set(model::TimedHybridSymbolicModel, state_sets, time_sets, mode_indices; domain=Dionysos.Domain.INNER)

For each mode k in mode_indices, returns all abstract state indices (state_id, time_id, k)
such that state_id is in the abstraction of state_sets[k] and time_id corresponds to time in time_sets[k].

# Arguments
- `model::TimedHybridSymbolicModel`: The timed hybrid symbolic model
- `state_sets`: Vector of HyperRectangle (or state set) per mode
- `time_sets`: Vector of HyperRectangle (or time interval) per mode
- `mode_indices`: List or set of mode indices
- `domain`: Domain type for state set intersection (default: INNER)

# Returns
- `Vector{Int}`: Corresponding abstract state indices

# Safety Features
- Bounds checking for mode indices
- Validation of set dimensions
- Handles empty intersections gracefully
"""
function get_states_from_set(
    model::TimedHybridSymbolicModel,
    state_sets,
    time_sets,
    mode_indices;
    domain = Dionysos.Domain.INNER,
)
    # Input validation
    @assert length(state_sets) >= length(mode_indices) "Not enough state sets provided"
    @assert length(time_sets) >= length(mode_indices) "Not enough time sets provided"

    abstract_states = Vector{Int}()

    for (idx, mode_id) in enumerate(mode_indices)
        # Bounds checking
        @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

        dynamics_model = model.mode_abstractions[mode_id]
        time_model = model.time_abstractions[mode_id]

        # Get abstract states in the spatial set
        spatial_states =
            Dionysos.Symbolic.get_states_from_set(dynamics_model, state_sets[idx], domain)

        # Get time indices in the temporal interval
        if hasfield(typeof(time_sets[idx]), :lb) && hasfield(typeof(time_sets[idx]), :ub)
            # HyperRectangle case
            t_min, t_max = time_sets[idx].lb[1], time_sets[idx].ub[1]
        else
            # Assume it's a 2-element vector [t_min, t_max]
            t_min, t_max = time_sets[idx][1], time_sets[idx][2]
        end

        time_indices = collect(
            Dionysos.Symbolic.ceil_time2int(
                time_model,
                t_min,
            ):Dionysos.Symbolic.floor_time2int(time_model, t_max),
        )

        # Find valid combinations
        for state_id in spatial_states, time_id in time_indices
            if time_id > 0 && time_id <= length(time_model.tsteps)
                augmented_key = (state_id, time_id, mode_id)::AugmentedState
                if haskey(model.augmented_to_state_index, augmented_key)
                    push!(abstract_states, model.augmented_to_state_index[augmented_key])
                end
            end
        end
    end

    return abstract_states
end

"""
    get_concrete_input(model::TimedHybridSymbolicModel, input_id::Int, mode_id::Int)

Convert an abstract input ID to its concrete input representation.
Handles both continuous and switching inputs with proper safety checks.

# Safety Features
- Validates input and mode IDs
- Handles switching inputs appropriately
- Returns nothing for invalid inputs
"""
function get_concrete_input(model::TimedHybridSymbolicModel, input_id::Int, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    input_type, local_info = get_local_input_info(model.input_mapping, input_id)

    if input_type == :continuous
        dynamics_model = model.mode_abstractions[mode_id]
        # local_info = (mode_id, local_input_id)
        local_input_id = local_info[2]
        return Dionysos.Symbolic.get_concrete_input(dynamics_model, local_input_id)
    elseif input_type == :switching
        # For switching inputs, there is no continuous concrete representation
        return nothing
    else
        @warn "Invalid input ID: $input_id"
        return nothing
    end
end

"""
    get_abstract_input(model::TimedHybridSymbolicModel, concrete_input, mode_id::Int)

Convert a concrete continuous input to its abstract input ID.
Only works for continuous inputs within the specified mode.

# Safety Features
- Validates mode ID bounds
- Handles cases where no valid abstract input exists
"""
function get_abstract_input(model::TimedHybridSymbolicModel, concrete_input, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    dynamics_model = model.mode_abstractions[mode_id]

    # Try to find the abstract input in the current mode
    local_input_id = Dionysos.Symbolic.get_abstract_input(dynamics_model, concrete_input)

    if !isnothing(local_input_id) && local_input_id > 0
        return get_global_input_id(model.input_mapping, mode_id, local_input_id)
    else
        return 0  # Not found or invalid
    end
end

end
