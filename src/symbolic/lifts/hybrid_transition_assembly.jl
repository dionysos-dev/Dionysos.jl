# ================================================================
# Transition assembly (pure: operates on already-built per-mode models)
# ================================================================

"""
    build_all_transitions(hs, mode_models, input_mapping) -> Vector{HybridTransition}

Assemble every hybrid transition: intra-mode (each mode model's own transitions,
which already encode any time advance) followed by inter-mode (guarded switches
between modes).
"""
function build_all_transitions(hs::HybridSystem, mode_models, input_mapping::GlobalInputMap)
    transition_list = Vector{HybridTransition}()
    sizehint!(transition_list, sum(get_n_transitions(m) for m in mode_models))

    add_intra_mode_transitions!(transition_list, mode_models, input_mapping)
    add_inter_mode_transitions!(transition_list, hs, mode_models, input_mapping)

    return transition_list
end

"Embed each mode model's own transitions into the global list, relabelling inputs."
function add_intra_mode_transitions!(
    transition_list,
    mode_models,
    input_mapping::GlobalInputMap,
)
    @assert !isempty(mode_models) "No mode models provided"

    for (mode_id, mode_model) in enumerate(mode_models)
        for (target_local, source_local, local_input_id) in enum_transitions(mode_model)
            global_input_id = get_global_input_id(input_mapping, mode_id, local_input_id)
            global_input_id > 0 || continue
            target_state = (target_local, mode_id)::HybridState
            source_state = (source_local, mode_id)::HybridState
            push!(
                transition_list,
                (target_state, source_state, global_input_id)::HybridTransition,
            )
        end
    end
end

"Add the guarded mode-switch transitions using each transition's guard and reset map."
function add_inter_mode_transitions!(
    transition_list,
    hs::HybridSystem,
    mode_models,
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

        @assert source_mode <= length(mode_models) "Source mode $source_mode out of bounds"
        @assert target_mode <= length(mode_models) "Target mode $target_mode out of bounds"

        reset_map = HybridSystems.resetmap(hs, transition)
        guard = HybridSystems.guard(hs, transition)

        if isnothing(guard)
            @warn "No guard found for transition $transition_id, skipping"
            continue
        end

        source_model = mode_models[source_mode]
        target_model = mode_models[target_mode]

        # Source states intersecting the guard (a box over the mode's `[x; t]` coord).
        source_locals = get_states_from_set(source_model, guard, MP.INNER)

        if isempty(source_locals)
            @warn "Empty guard intersection for transition $transition_id, skipping"
            continue
        end

        for source_local in source_locals
            # Concretize, apply the reset map, and abstract in the target mode.
            source_coord = get_concrete_state(source_model, source_local)
            reset_result = MS.apply(reset_map, source_coord)
            target_local = abstract_switch_target(target_model, reset_result)

            target_local > 0 || continue
            target_state = (target_local, target_mode)::HybridState
            source_state = (source_local, source_mode)::HybridState
            push!(
                transition_list,
                (target_state, source_state, global_input_id)::HybridTransition,
            )
        end
    end
end

"""
    build_symbolic_automaton(transition_list, mode_models, input_mapping) -> (flat, automaton)

Flatten the hybrid transitions into an [`IndexedAutomatonList`](@ref), returning
`(flat, automaton)` where `flat` is the [`FlatIndex`](@ref) between the integer
numbering and the `(local_id, mode_id)` pairs.
"""
function build_symbolic_automaton(
    transition_list,
    mode_models,
    input_mapping::GlobalInputMap,
)
    @assert !isempty(transition_list) "Transition list cannot be empty"
    @assert !isempty(mode_models) "Mode models cannot be empty"

    estimated_states = sum(get_n_state(m) for m in mode_models)

    states = Set{HybridState}()
    inputs_set = Set{Int}()
    sizehint!(states, estimated_states)
    sizehint!(inputs_set, input_mapping.total_inputs)

    for (target, source, input) in transition_list
        push!(states, target)
        push!(states, source)
        push!(inputs_set, input)
    end

    flat = FlatIndex(collect(states))
    ninputs = length(inputs_set)

    symbolic_automaton = IndexedAutomatonList(n_flat(flat), ninputs)

    @inbounds for (target, source, abstract_input) in transition_list
        target_int = flat_id(flat, target)
        source_int = flat_id(flat, source)
        add_transition!(symbolic_automaton, source_int, target_int, abstract_input)
    end

    return flat, symbolic_automaton
end
