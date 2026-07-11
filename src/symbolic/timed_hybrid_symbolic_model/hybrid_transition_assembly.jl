# ================================================================
# Guard geometry helpers (guard assumed to be a box: [spatial…, time])
# ================================================================

"Spatial part (all but the last dimension) of a box `guard`."
function extract_spatial_part(guard)
    if isa(guard, LazySets.AbstractHyperrectangle)
        return UT.box(LazySets.low(guard)[1:(end - 1)], LazySets.high(guard)[1:(end - 1)])
    else
        error("Unsupported guard type: $(typeof(guard))")
    end
end

"Temporal part `[t_min, t_max]` (last dimension) of a box `guard`."
function extract_temporal_part(guard)
    if isa(guard, LazySets.AbstractHyperrectangle)
        return [LazySets.low(guard)[end], LazySets.high(guard)[end]]
    else
        error("Unsupported guard type: $(typeof(guard))")
    end
end

"Time indices of `time_model` falling within the `[t_min, t_max]` interval."
function get_time_indices_from_interval(time_model, temporal_interval)
    t_min, t_max = temporal_interval
    return findall(t -> t_min <= t <= t_max, time_model.tsteps)
end

# ================================================================
# Transition assembly (pure: operates on already-built mode abstractions)
# ================================================================

"""
    build_all_transitions(hs, mode_abstractions, input_mapping) -> Vector{TransitionTuple}

Assemble every augmented transition of the timed hybrid system: intra-mode
(dynamics × time steps within each mode) followed by inter-mode (guarded
switches between modes).
"""
function build_all_transitions(
    hs::HybridSystem,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)
    intra_mode_transitions = sum(
        get_n_transitions(abs_pair[1]) * length(abs_pair[2].tsteps) for
        abs_pair in mode_abstractions
    )

    transition_list = Vector{TransitionTuple}()
    sizehint!(transition_list, intra_mode_transitions)

    add_intra_mode_transitions!(transition_list, mode_abstractions, input_mapping)
    add_inter_mode_transitions!(transition_list, hs, mode_abstractions, input_mapping)

    return transition_list
end

"Add the within-mode transitions (spatial dynamics advanced across each time step)."
function add_intra_mode_transitions!(
    transition_list,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)
    n_modes = length(mode_abstractions)
    @assert n_modes > 0 "No mode abstractions provided"

    for (mode_id, (symmodel_dynamics, symmodel_time)) in enumerate(mode_abstractions)
        time_steps = symmodel_time.tsteps
        n_time_steps = length(time_steps)

        if n_time_steps == 1
            @inbounds for (target, source, local_input_id) in
                          enum_transitions(symmodel_dynamics)
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
            transitions_cache = collect(enum_transitions(symmodel_dynamics))
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

"Add the guarded mode-switch transitions using each transition's guard and reset map."
function add_inter_mode_transitions!(
    transition_list,
    hs::HybridSystem,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)
    transitions = collect(HybridSystems.transitions(hs.automaton))

    for (transition_id, transition) in enumerate(transitions)
        global_input_id = get_switching_global_id(input_mapping, transition_id)

        if global_input_id <= 0
            @warn "Invalid global input ID for transition $transition_id, skipping"
            continue
        end

        source_mode = HybridSystems.source(hs.automaton, transition)
        target_mode = HybridSystems.target(hs.automaton, transition)

        @assert source_mode <= length(mode_abstractions) "Source mode $source_mode out of bounds"
        @assert target_mode <= length(mode_abstractions) "Target mode $target_mode out of bounds"

        reset_map = HybridSystems.resetmap(hs, transition)
        guard = HybridSystems.guard(hs, transition)

        if isnothing(guard)
            @warn "No guard found for transition $transition_id, skipping"
            continue
        end

        (source_symmodel_dynamics, source_time_model) = mode_abstractions[source_mode]
        (target_symmodel_dynamics, target_time_model) = mode_abstractions[target_mode]

        # Split the guard into spatial and temporal parts
        guard_spatial = extract_spatial_part(guard)
        guard_temporal = extract_temporal_part(guard)

        # Source states / times intersecting the guard
        source_states =
            get_states_from_set(source_symmodel_dynamics, guard_spatial, MP.INNER)
        time_indices = get_time_indices_from_interval(source_time_model, guard_temporal)

        if isempty(source_states) || isempty(time_indices)
            @warn "Empty guard intersection for transition $transition_id, skipping"
            continue
        end

        for source_state in source_states, source_time_idx in time_indices
            if source_time_idx > length(source_time_model.tsteps) || source_time_idx <= 0
                continue
            end

            # Build the augmented source state [x1, …, xn, t], apply the reset map
            source_continuous_state =
                get_concrete_state(source_symmodel_dynamics, source_state)
            source_time_value = source_time_model.tsteps[source_time_idx]
            augmented_source_state = vcat(source_continuous_state, source_time_value)

            reset_result = MS.apply(reset_map, augmented_source_state)
            reset_continuous_part = reset_result[1:(end - 1)]
            reset_time_value = reset_result[end]

            # Corresponding target symbolic state and time index
            target_state =
                find_symbolic_state(target_symmodel_dynamics, reset_continuous_part)
            target_time_idx = ceil_time2int(target_time_model, reset_time_value)

            if target_state > 0 &&
               target_time_idx > 0 &&
               target_time_idx <= length(target_time_model.tsteps)
                target_aug_state =
                    (target_state, target_time_idx, target_mode)::AugmentedState
                source_aug_state =
                    (source_state, source_time_idx, source_mode)::AugmentedState
                push!(
                    transition_list,
                    (target_aug_state, source_aug_state, global_input_id)::TransitionTuple,
                )
            end
        end
    end
end

"""
    build_symbolic_automaton(transition_list, mode_abstractions, input_mapping)

Flatten the augmented transitions into an [`IndexedAutomatonList`](@ref),
returning `(state_index_to_augmented, augmented_to_state_index, automaton)`.
"""
function build_symbolic_automaton(
    transition_list,
    mode_abstractions,
    input_mapping::GlobalInputMap,
)
    @assert !isempty(transition_list) "Transition list cannot be empty"
    @assert !isempty(mode_abstractions) "Mode abstractions cannot be empty"

    estimated_states = sum(
        get_n_state(abs_pair[1]) * length(abs_pair[2].tsteps) for
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

    symbolic_automaton = IndexedAutomatonList(nstates, ninputs)

    @inbounds for (target, source, abstract_input) in transition_list
        target_int = augmented_to_state_index[target]
        source_int = augmented_to_state_index[source]

        add_transition!(symbolic_automaton, source_int, target_int, abstract_input)
    end

    return state_index_to_augmented, augmented_to_state_index, symbolic_automaton
end
