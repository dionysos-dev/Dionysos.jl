# A hybrid state pairs a per-mode local state id with a mode id: `(local_id, mode_id)`.
# The local id is opaque to the composition: for a clock-lifted mode it already
# encodes `(spatial_state, time_id)`; for a plain mode it is just the spatial state.
const HybridState = Tuple{Int, Int}
const HybridTransition = Tuple{HybridState, HybridState, Int}

"""
    HybridSymbolicModel{Mods, A, G} <: AbstractSymbolicModel

Symbolic abstraction of a hybrid system, composed from one `AbstractSymbolicModel`
per mode via a [`ModeLift`](@ref). Each mode model is either a plain spatial
abstraction (`x`) or a [`ClockLiftedSymbolicModel`](@ref) (`(x, t)`); the
`(local_state_id, mode_id)` pairs are flattened to a single integer numbering and
wired into one automaton, with inputs unified through a [`GlobalInputMap`](@ref).

# Fields
- `mode_models`: per-mode symbolic models.
- `flat`: [`FlatIndex`](@ref) bijection between the integer numbering and the
  `(local_state_id, mode_id)` pairs.
- `automaton`: the flattened transition automaton.
- `input_mapping`: the global input map.
"""
struct HybridSymbolicModel{Mods, A, G} <: AbstractSymbolicModel
    mode_models::Mods
    flat::FlatIndex{HybridState}
    automaton::A
    input_mapping::G
end

# Transitional alias while the hybrid optimizers/tests migrate off the old name.
const TimedHybridSymbolicModel = HybridSymbolicModel

get_automaton(sym::HybridSymbolicModel) = sym.automaton

# ================================================================
# Accessors
# ================================================================

"Number of (flattened) hybrid states."
get_n_state(model::HybridSymbolicModel) = n_flat(model.flat)

"Total number of global inputs (continuous + switching)."
get_n_input(model::HybridSymbolicModel) = model.input_mapping.total_inputs

enum_states(model::HybridSymbolicModel) = 1:get_n_state(model)

"Enumerate the local input ids of `mode_id`."
function enum_inputs(model::HybridSymbolicModel, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_models) "Mode ID $mode_id out of bounds"
    return enum_inputs(model.mode_models[mode_id])
end

"""
    get_concrete_state(model::HybridSymbolicModel, state_index) -> (x, t, mode_id)

Concretize a flattened state index to `(continuous_state, time_value, mode_id)`.
"""
function get_concrete_state(model::HybridSymbolicModel, state_index::Int)
    @assert 1 <= state_index <= n_flat(model.flat) "State index $state_index out of bounds"
    (local_id, mode_id) = flat_key(model.flat, state_index)
    @assert 1 <= mode_id <= length(model.mode_models) "Mode ID $mode_id out of bounds"
    x, t = base_state_and_time(model.mode_models[mode_id], local_id)
    return (x, t, mode_id)
end

"""
    get_abstract_state(model::HybridSymbolicModel, (x, t, mode_id)) -> Int

Abstract a concrete augmented state to its flattened index (`0` if the augmented
key is not present in the model).
"""
function get_abstract_state(model::HybridSymbolicModel, augmented_state)
    (x, t, mode_id) = augmented_state
    @assert 1 <= mode_id <= length(model.mode_models) "Mode ID $mode_id out of bounds"
    local_id = get_abstract_state(model.mode_models[mode_id], vcat(x, SVector(t)))
    local_id <= 0 &&
        error("No valid abstract state found for augmented_state $augmented_state")
    return flat_id(model.flat, (local_id, mode_id))
end

"""
    get_states_from_set(model, state_sets, time_sets, mode_indices; domain=MP.INNER)

For each mode in `mode_indices`, collect the flattened state indices whose spatial
part lies in `state_sets[idx]` and whose time index lies in `time_sets[idx]`.
"""
function get_states_from_set(
    model::HybridSymbolicModel,
    state_sets,
    time_sets,
    mode_indices;
    domain = MP.INNER,
)
    @assert length(state_sets) >= length(mode_indices) "Not enough state sets provided"
    @assert length(time_sets) >= length(mode_indices) "Not enough time sets provided"

    abstract_states = Vector{Int}()
    for (idx, mode_id) in enumerate(mode_indices)
        @assert 1 <= mode_id <= length(model.mode_models) "Mode ID $mode_id out of bounds"
        combined = _augment_box(state_sets[idx], time_sets[idx])
        local_ids = get_states_from_set(model.mode_models[mode_id], combined, domain)
        for local_id in local_ids
            gid = flat_id(model.flat, (local_id, mode_id))
            gid > 0 && push!(abstract_states, gid)
        end
    end
    return abstract_states
end

# Combine a spatial box and a time box/interval into a single `[x; t]` box.
function _augment_box(state_set, time_set)
    if isa(time_set, LazySets.AbstractHyperrectangle)
        tmin = LazySets.low(time_set, 1)
        tmax = LazySets.high(time_set, 1)
    else
        tmin, tmax = time_set[1], time_set[2]
    end
    return UT.box(vcat(LazySets.low(state_set), tmin), vcat(LazySets.high(state_set), tmax))
end

"Concretize global input `input_id` in `mode_id` (`nothing` for switching inputs)."
function get_concrete_input(model::HybridSymbolicModel, input_id::Int, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_models) "Mode ID $mode_id out of bounds"

    input_type, local_info = get_local_input_info(model.input_mapping, input_id)

    if input_type == :continuous
        local_input_id = local_info[2]
        return get_concrete_input(model.mode_models[mode_id], local_input_id)
    elseif input_type == :switching
        return nothing
    else
        @warn "Invalid input ID: $input_id"
        return nothing
    end
end

"Abstract a concrete input in `mode_id` to its global input id (`0` if not found)."
function get_abstract_input(model::HybridSymbolicModel, concrete_input, mode_id::Int)
    @assert 1 <= mode_id <= length(model.mode_models) "Mode ID $mode_id out of bounds"

    local_input_id = get_abstract_input(model.mode_models[mode_id], concrete_input)

    if !isnothing(local_input_id) && local_input_id > 0
        return get_global_input_id(model.input_mapping, mode_id, local_input_id)
    else
        return 0
    end
end
