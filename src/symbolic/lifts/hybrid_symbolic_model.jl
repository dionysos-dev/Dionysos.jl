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

get_automaton(sym::HybridSymbolicModel) = sym.automaton

"The [`GlobalInputMap`](@ref) unifying the modes' inputs and the mode switches."
get_global_input_map(sym::HybridSymbolicModel) = sym.input_mapping

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

# A mode model contributes its own concrete coordinate to the hybrid augmented
# state: `(x, t)` when clock-lifted, `(x,)` when time-free. The hybrid appends the
# mode id, so the augmented state is `(x, t, mode)` or `(x, mode)` accordingly.
mode_concrete_coord(m::ClockLiftedSymbolicModel, local_id::Int) =
    base_state_and_time(m, local_id)
mode_concrete_coord(m::SymbolicModel, local_id::Int) = (get_concrete_state(m, local_id),)

# Local abstract state id of a mode model from the hybrid augmented state (0 if absent).
_local_abstract_state(m::ClockLiftedSymbolicModel, aug) =
    get_abstract_state(m, vcat(aug[1], SVector(aug[2])))
function _local_abstract_state(m::SymbolicModel, aug)
    q = get_abstract_state(m, aug[1])
    return (q === nothing || q <= 0) ? 0 : q
end

"""
    get_concrete_state(model::HybridSymbolicModel, state_index) -> (x[, t], mode_id)

Concretize a flattened state index to the augmented concrete state: `(x, t, mode)`
for a clock-lifted mode, `(x, mode)` for a time-free mode.
"""
function get_concrete_state(model::HybridSymbolicModel, state_index::Int)
    @assert 1 <= state_index <= n_flat(model.flat) "State index $state_index out of bounds"
    (local_id, mode_id) = flat_key(model.flat, state_index)
    @assert 1 <= mode_id <= length(model.mode_models) "Mode ID $mode_id out of bounds"
    return (mode_concrete_coord(model.mode_models[mode_id], local_id)..., mode_id)
end

"""
    get_abstract_state(model::HybridSymbolicModel, augmented_state) -> Int

Abstract a concrete augmented state (`(x, t, mode)` or `(x, mode)`) to its flattened
index, or `0` if it is out of the abstraction's domain (unknown mode, or a
spatial/time coordinate outside the grid). Callers that require a valid state
(e.g. an initial state) should check for `0`.
"""
function get_abstract_state(model::HybridSymbolicModel, augmented_state)
    mode_id = augmented_state[end]
    (1 <= mode_id <= length(model.mode_models)) || return 0
    local_id = _local_abstract_state(model.mode_models[mode_id], augmented_state)
    local_id <= 0 && return 0
    return flat_id(model.flat, (local_id, mode_id))
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
