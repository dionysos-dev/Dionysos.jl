# An augmented state pairs a spatial abstract state with a time index and a mode:
# `(state_id, time_id, mode_id)`. Transitions carry a global input id.
const AugmentedState = Tuple{Int, Int, Int}
const TransitionTuple = Tuple{AugmentedState, AugmentedState, Int}

"""
    TimedHybridSymbolicModel{S1, A, T, G} <: AbstractSymbolicModel

Symbolic abstraction of a timed hybrid system. Each mode contributes a spatial
dynamics abstraction and a [`TimeSymbolicModel`](@ref); the augmented states
`(state_id, time_id, mode_id)` are flattened to a single integer numbering and
wired into one automaton, with inputs unified through a [`GlobalInputMap`](@ref).

# Fields
- `mode_abstractions`: per-mode spatial dynamics symbolic models.
- `time_abstractions`: per-mode time symbolic models.
- `state_index_to_augmented` / `augmented_to_state_index`: integer ↔ augmented state.
- `symbolic_automaton`: the flattened transition automaton.
- `input_mapping`: the global input map.
"""
struct TimedHybridSymbolicModel{S1, A, T, G} <: AbstractSymbolicModel
    mode_abstractions::Vector{S1}
    time_abstractions::Vector{T}
    state_index_to_augmented::Vector{AugmentedState}
    augmented_to_state_index::Dict{AugmentedState, Int}
    symbolic_automaton::A
    input_mapping::G
end

get_automaton(sym::TimedHybridSymbolicModel) = sym.symbolic_automaton

# ================================================================
# Accessors
# ================================================================

"Number of (flattened) augmented states."
get_n_state(model::TimedHybridSymbolicModel) = length(model.state_index_to_augmented)

"Total number of global inputs (continuous + switching)."
get_n_input(model::TimedHybridSymbolicModel) = model.input_mapping.total_inputs

enum_states(model::TimedHybridSymbolicModel) = 1:get_n_state(model)

"Enumerate the local input ids of `mode_id`."
function enum_inputs(model::TimedHybridSymbolicModel, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"
    return enum_inputs(model.mode_abstractions[mode_id])
end

"""
    get_concrete_state(model::TimedHybridSymbolicModel, state_index) -> (x, t, mode_id)

Concretize a flattened state index to `(continuous_state, time_value, mode_id)`.
"""
function get_concrete_state(model::TimedHybridSymbolicModel, state_index::Int)
    @assert 1 <= state_index <= length(model.state_index_to_augmented) "State index $state_index out of bounds"

    (state_id, time_id, mode_id) = model.state_index_to_augmented[state_index]

    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    dynamics_model = model.mode_abstractions[mode_id]
    time_model = model.time_abstractions[mode_id]

    continuous_state = get_concrete_state(dynamics_model, state_id)
    time_value = int2time(time_model, time_id)

    return (continuous_state, time_value, mode_id)
end

"""
    get_abstract_state(model::TimedHybridSymbolicModel, (x, t, mode_id)) -> Int

Abstract a concrete augmented state to its flattened index (`0` if the augmented
key is not present in the model).
"""
function get_abstract_state(model::TimedHybridSymbolicModel, augmented_state)
    (continuous_state, time_value, mode_id) = augmented_state

    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    dynamics_model = model.mode_abstractions[mode_id]
    time_model = model.time_abstractions[mode_id]

    state_id = get_abstract_state(dynamics_model, continuous_state)
    time_id = floor_time2int(time_model, time_value)

    if isnothing(state_id) || time_id <= 0
        error("No valid abstract state found for augmented_state $augmented_state")
    end

    augmented_key = (state_id, time_id, mode_id)::AugmentedState
    return get(model.augmented_to_state_index, augmented_key, 0)
end

"""
    get_states_from_set(model, state_sets, time_sets, mode_indices; domain=MP.INNER)

For each mode in `mode_indices`, collect the flattened state indices whose spatial
part lies in `state_sets[idx]` and whose time index lies in `time_sets[idx]`.
"""
function get_states_from_set(
    model::TimedHybridSymbolicModel,
    state_sets,
    time_sets,
    mode_indices;
    domain = MP.INNER,
)
    @assert length(state_sets) >= length(mode_indices) "Not enough state sets provided"
    @assert length(time_sets) >= length(mode_indices) "Not enough time sets provided"

    abstract_states = Vector{Int}()

    for (idx, mode_id) in enumerate(mode_indices)
        @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

        dynamics_model = model.mode_abstractions[mode_id]
        time_model = model.time_abstractions[mode_id]

        spatial_states = get_states_from_set(dynamics_model, state_sets[idx], domain)

        if isa(time_sets[idx], LazySets.AbstractHyperrectangle)
            t_min = LazySets.low(time_sets[idx], 1)
            t_max = LazySets.high(time_sets[idx], 1)
        else
            t_min, t_max = time_sets[idx][1], time_sets[idx][2]
        end

        time_indices =
            collect(ceil_time2int(time_model, t_min):floor_time2int(time_model, t_max))

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

"Concretize global input `input_id` in `mode_id` (`nothing` for switching inputs)."
function get_concrete_input(model::TimedHybridSymbolicModel, input_id::Int, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    input_type, local_info = get_local_input_info(model.input_mapping, input_id)

    if input_type == :continuous
        dynamics_model = model.mode_abstractions[mode_id]
        local_input_id = local_info[2]
        return get_concrete_input(dynamics_model, local_input_id)
    elseif input_type == :switching
        return nothing
    else
        @warn "Invalid input ID: $input_id"
        return nothing
    end
end

"Abstract a concrete input in `mode_id` to its global input id (`0` if not found)."
function get_abstract_input(model::TimedHybridSymbolicModel, concrete_input, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_abstractions) "Mode ID $mode_id out of bounds"

    dynamics_model = model.mode_abstractions[mode_id]
    local_input_id = get_abstract_input(dynamics_model, concrete_input)

    if !isnothing(local_input_id) && local_input_id > 0
        return get_global_input_id(model.input_mapping, mode_id, local_input_id)
    else
        return 0
    end
end

"""
    find_symbolic_state(symmodel, continuous_state) -> Int

Abstract-state index of `continuous_state` in a per-mode dynamics model, or `0`
if the state is empty/`nothing` or has no valid abstraction.
"""
function find_symbolic_state(symmodel, continuous_state)
    if isnothing(continuous_state) || isempty(continuous_state)
        return 0
    end

    state_idx = get_abstract_state(symmodel, continuous_state)
    return (isnothing(state_idx) || state_idx <= 0) ? 0 : state_idx
end
