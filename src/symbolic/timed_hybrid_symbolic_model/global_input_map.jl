"""
    GlobalInputMap

Bijection between per-mode local input ids and a single global input alphabet
for a timed hybrid system. Continuous inputs (one block per mode) occupy
`continuous_range`; switching inputs (one per hybrid-automaton transition)
occupy `switching_range`.

# Fields
- `total_inputs`, `continuous_inputs`, `switching_inputs`: counts.
- `continuous_to_global` / `global_to_continuous`: `(mode_id, local_input_id)` ↔ `global_id`.
- `switching_to_global` / `global_to_switching`: `transition_id` ↔ `global_id`.
- `continuous_range`, `switching_range`: the id ranges of each block.
- `switch_labels`: human-readable `"SWITCH src -> tgt"` label per switching input.
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
    switch_labels::Vector{String}
end

"""
    GlobalInputMap(mode_abstractions, hs::HybridSystem)

Build the [`GlobalInputMap`](@ref) for a hybrid system from its per-mode
`(symbolic_dynamics, symbolic_time)` abstractions: continuous inputs are laid
out mode by mode, then one switching input per hybrid-automaton transition.
"""
function GlobalInputMap(abstract_systems, hs::HybridSystem)
    # Step 1: Allocate continuous inputs
    continuous_to_global = Dict{Tuple{Int, Int}, Int}()
    global_to_continuous = Dict{Int, Tuple{Int, Int}}()
    continuous_count = 0
    for (mode_id, (symmodel_dynam, _)) in enumerate(abstract_systems)
        input_count = get_n_input(symmodel_dynam)
        for local_input_id in 1:input_count
            global_id = continuous_count + local_input_id
            continuous_to_global[(mode_id, local_input_id)] = global_id
            global_to_continuous[global_id] = (mode_id, local_input_id)
        end
        continuous_count += input_count
    end
    # Step 2: Allocate switching inputs and labels
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
    # Step 3: Compute ranges
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

# === Accessors ===

"Global input id of a local continuous input `(mode_id, local_input_id)` (`0` if absent)."
function get_global_input_id(gim::GlobalInputMap, mode_id::Int, local_input_id::Int)
    return get(gim.continuous_to_global, (mode_id, local_input_id), 0)
end

"Global input id of the switching input for `transition_id` (`0` if absent)."
function get_switching_global_id(gim::GlobalInputMap, transition_id::Int)
    return get(gim.switching_to_global, transition_id, 0)
end

"""
    get_local_input_info(gim::GlobalInputMap, global_id) -> (kind, info)

Classify a global input id: returns `(:continuous, (mode_id, local_input_id))`,
`(:switching, transition_id)`, or `(:invalid, nothing)`.
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

"Whether `global_id` is a continuous input."
function is_continuous_input(gim::GlobalInputMap, global_id::Int)
    return global_id in gim.continuous_range
end

"Whether `global_id` is a switching input."
function is_switching_input(gim::GlobalInputMap, global_id::Int)
    return global_id in gim.switching_range
end
