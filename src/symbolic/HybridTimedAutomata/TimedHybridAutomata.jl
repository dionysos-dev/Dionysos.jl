export TimedHybridAutomata

module TimedHybridAutomata

using HybridSystems, MathematicalSystems, Dionysos, StaticArrays

# actuellement un problème de dépendance circulaire -> LoadError: UndefVarError: `Problem` not defined in `Dionysos`
# Actuellement j'utilise la fonction `Build_Timed_Hybrid_Automaton` pour construire le modèle symbolique hybride temporel complet.
# A terme, on devrait permettre à l'utilisateur de fournir ses propres méthodes d'abstraction et les optimizers qu'il souhaite utiliser pour chaque mode.
# Cependant, cela créé un problème de dépendance circulaire entre les modules `HybridTimedAutomata`, 'optim' et `Problem`.
import Dionysos
#PR = Dionysos.Problem
#OP = Dionysos.Optim
#AB = OP.Abstraction


# ================================================================
# Symbolic temporal hybrid model structure
# ================================================================

"""
    TemporalHybridSymbolicModel{S1, A, N, T, G}

Main structure representing the symbolic abstraction of a temporal hybrid system.

# Fields
- `symmodels::Vector{S1}`: Symbolic models of the dynamics per mode
- `time_symbolic_models::Vector{T}`: Symbolic models of time per mode
- `int2aug_state::Vector{NTuple{N, Int}}`: Integer → augmented state (state_symbol, time_symbol, mode_id)
- `aug_state2int::Dict{NTuple{N, Int}, Int}`: Augmented state → integer
- `autom::A`: Final symbolic automaton
- `global_input_map::G`: Unified global input system
"""
mutable struct TemporalHybridSymbolicModel{S1, A, N, T, G}
    symmodels::Vector{S1}
    time_symbolic_models::Vector{T}
    int2aug_state::Vector{NTuple{N, Int}}
    aug_state2int::Dict{NTuple{N, Int}, Int}
    autom::A
    global_input_map::G
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
        input_count = SY.get_n_input(symmodel_dynam)
        for local_input_id in 1:input_count
            global_id = continuous_count + local_input_id
            continuous_to_global[(mode_id, local_input_id)] = global_id
            global_to_continuous[global_id] = (mode_id, local_input_id)
        end
        continuous_count += input_count
    end
    # Phase 2: Allocate switching inputs
    switching_to_global = Dict{Int, Int}()
    global_to_switching = Dict{Int, Int}()
    switching_count = 0
    for (transition_id, _) in enumerate(HybridSystems.transitions(hs.automaton))
        global_id = continuous_count + switching_count + 1
        switching_to_global[transition_id] = global_id
        global_to_switching[global_id] = transition_id
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
        switching_range
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

# Arguments
- `system`: The continuous system to abstract (currently mainly ConstrainedBlackBoxControlContinuousSystem)
- `growth_bound`: Bound on the system's growth (for over-approximation)
- `param_discretisation`: Tuple (dx, du, tstep) for state, input, and time discretization

# Returns
- Symbolic abstraction of the system
"""
function build_dynamical_symbolic_model(system, growth_bound, param_discretisation)
    
    # Build the symbolic model for the dynamics
    problem = PR.EmptyProblem(system, system.X)
    dx, du, tstep = param_discretisation
    nx = system.statedim
    nu = system.inputdim
    x0 = SVector{nx, Float64}(zeros(nx))
    hx = SVector{nx, Float64}(fill(dx, nx))
    state_grid = DO.GridFree(x0, hx)
    u0 = SVector{nu, Float64}(zeros(nu))
    hu = SVector{nu, Float64}(fill(du, nu))
    input_grid = DO.GridFree(u0, hu)
    function my_jacobian_bound(u)
        return growth_bound
    end
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("approx_mode"), AB.UniformGridAbstraction.GROWTH)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), my_jacobian_bound)
    MOI.optimize!(optimizer)
    return MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
end

"""
    build_initial_symmodel_by_mode(hs::HybridSystem, growth_bounds::SVector{}, param_discretisation)

Build a list of symbolic models (dynamics and time) for each mode of a hybrid system.

# Arguments
- `hs::HybridSystem`: The hybrid system
- `growth_bounds::SVector{}`: Growth bounds per mode
- `param_discretisation`: Discretization parameters per mode

# Returns
- Vector of (symbolic_dynamics, symbolic_time) tuples per mode
"""
function build_initial_symmodel_by_mode(hs::HybridSystem, growth_bounds::SVector{}, param_discretisation)
    # Build a list of symbolic models for each mode of the hybrid system
    abstract_systems = []
    for (i, mode_id) in enumerate(HybridSystems.states(hs.automaton))
        mode_system = HybridSystems.mode(hs, mode_id)
        dyn_sys = mode_system.systems[1]    # physical dynamics
        time_sys = mode_system.systems[2]   # time dynamics
        symmodel_dynam = build_dynamical_symbolic_model(dyn_sys, growth_bounds[i], param_discretisation[i])
        symmodel_time = BuildTimeSymbolicModel(time_sys, param_discretisation[i][3])
        push!(abstract_systems, (symmodel_dynam, symmodel_time))
    end
    return abstract_systems
end


# ================================================================
# Functions to add transitions
# ================================================================

"""
    add_mode_transitions!(translist, abstract_systems, input_mapping::GlobalInputMap)

Add intra-mode transitions to the transition list, based on the symbolic models and time discretization.

# Arguments
- `translist`: The list to which transitions are appended
- `abstract_systems`: Vector of (symbolic_dynamics, symbolic_time) tuples per mode
- `input_mapping::GlobalInputMap`: The global input mapping
"""
function add_mode_transitions!(translist, abstract_systems, input_mapping::GlobalInputMap)
    # Add transitions for each mode based on the symbolic models (cartesian product of states and time in the case of time taken into account)
    for (mode_id, (symmodel_dynam, symmodel_time)) in enumerate(abstract_systems)
        tm = symmodel_time
        abstract_system = symmodel_dynam
        for (target, source, local_input_id) in SY.enum_transitions(abstract_system)
            global_input_id = get_global_input_id(input_mapping, mode_id, local_input_id)
            if length(tm.tsteps) == 1
                # Time-frozen: only one transition (t=1 symbolically)
                push!(translist, ((target, 1, mode_id), (source, 1, mode_id), global_input_id))
            else
                for k in 1:length(tm.tsteps) - 1
                    # For each time step, add a transition ((target, k+1, mode_id), (source, k, mode_id), global_input_id) with k being the time index of the source state
                    push!(translist, ((target, k + 1, mode_id), (source, k, mode_id), global_input_id))
                end
            end
        end
    end
end

# Note
# This function is "not robust" for all possible reset maps. If the new state value obtained after applying the reset map lies exactly on the boundary between two cells in the target mode's state space, the algorithm may fail to find a valid symbolic state.  
# It is necessary to ensure that the reset map correctly maps all values from the guard of the source mode to existing values in the state space of the target mode.  
# I will improve this function in the future.

"""
    add_switching_transitions!(translist, hs::HybridSystem, abstract_systems, input_mapping::GlobalInputMap)

Add switching transitions (between modes) to the transition list, using guards and reset maps.

# Arguments
- `translist`: The list to which transitions are appended
- `hs::HybridSystem`: The hybrid system
- `abstract_systems`: Vector of (symbolic_dynamics, symbolic_time) tuples per mode
- `input_mapping::GlobalInputMap`: The global input mapping
"""
function add_switching_transitions!(translist, hs::HybridSystem, abstract_systems, input_mapping::GlobalInputMap)
    for (transition_id, transition) in enumerate(HybridSystems.transitions(hs.automaton))
        # Get the global input id for this switching transition
        global_input_id = get_switching_global_id(input_mapping, transition_id)

        # Extract source and target mode indices
        source_mode = HybridSystems.source(hs.automaton, transition)
        target_mode = HybridSystems.target(hs.automaton, transition)

        # Get the reset map and guard for this transition
        reset_map = HybridSystems.resetmap(hs, transition)
        guard = HybridSystems.guard(hs, transition)

        # Get the symbolic models for the source and target modes
        (source_symmodel_dynam, source_tm) = abstract_systems[source_mode]
        (target_symmodel_dynam, target_tm) = abstract_systems[target_mode]

        # Split the guard into spatial and temporal parts
        guard_spatial = extract_spatial_part(guard)
        guard_temporal = extract_temporal_part(guard)

        # Get all source states that intersect with the spatial guard
        source_states = SY.get_states_from_set(source_symmodel_dynam, guard_spatial, DO.INNER)
        # Get all time indices that intersect with the temporal guard
        time_indices = get_time_indices_from_interval(source_tm, guard_temporal)

        # For each combination of (state, time) in the guard
        for source_state in source_states, source_time_idx in time_indices
            # Build the augmented source state [x1, x2, ..., xn, t]
            source_continuous_state = SY.get_concrete_state(source_symmodel_dynam, source_state)
            source_time_value = source_tm.tsteps[source_time_idx]
            augmented_source_state = vcat(source_continuous_state, source_time_value)

            # Apply the reset map to the augmented state
            reset_result = MathematicalSystems.apply(reset_map, augmented_source_state)
            reset_continuous_part = reset_result[1:end-1]
            reset_time_value = reset_result[end]

            # Find the corresponding target symbolic state and time index
            target_state = find_symbolic_state(target_symmodel_dynam, reset_continuous_part)
            target_time_idx = find_time_index(target_tm, reset_time_value)

            # Add the transition if both target state and time are valid
            if target_state > 0 && target_time_idx > 0
                push!(translist, (
                    (target_state, target_time_idx, target_mode),
                    (source_state, source_time_idx, source_mode),
                    global_input_id
                ))
            end
        end
    end
end

# ================================================================
# Build the final automaton from the transition list
# ================================================================

"""
    build_automaton(translist)

Build the symbolic automaton from the list of temporal transitions.

# Arguments
- `translist`: List of temporal transitions (tuples of (target, source, input))

# Returns
- Tuple (int2aug_state, aug_state2int, automaton):
    - `int2aug_state`: Vector mapping integer indices to augmented states
    - `aug_state2int`: Dict mapping augmented states to integer indices
    - `automaton`: The constructed symbolic automaton
"""
function build_automaton(translist)
    function enum_augmented_states(transitions)
        states = Set{Any}()
        for (target, source, _) in transitions
            push!(states, target)
            push!(states, source)
        end
        return collect(states)
    end
    augmented_states = enum_augmented_states(translist)
    nstates = length(augmented_states)
    int2aug_state = [aug_state for aug_state in augmented_states]
    aug_state2int = Dict((aug_state, i) for (i, aug_state) in enumerate(augmented_states))
    inputs_set = Set{Int}()
    for (_, _, input) in translist
        push!(inputs_set, input)
    end
    ninputs = length(inputs_set)
    autom = SY.NewIndexedAutomatonList(nstates, ninputs)
    for (target, source, abstract_input) in translist
        target_int = aug_state2int[target]
        source_int = aug_state2int[source]
        SY.add_transition!(autom, source_int, target_int, abstract_input)
    end
    return int2aug_state, aug_state2int, autom
end

# ================================================================
# Temporal hybrid symbolic model constructor
# ================================================================

# C'est la fonction que l'utilisateur doit appeler pour construire le modèle symbolique hybride temporel complet.
# A terme comme précisé ci dessus, on devrait permettre à l'utilisateur de fournir ses propres méthodes d'abstraction et les optimizers qu'il souhaite utiliser pour chaque mode.

"""
    Build_Timed_Hybrid_Automaton(hs::HybridSystem, growth_bounds::SVector{}, param_discretisation)

Construct the full temporal hybrid symbolic model for a given hybrid system.

# Arguments
- `hs::HybridSystem`: The hybrid system
- `growth_bounds::SVector{}`: Growth bounds per mode
- `param_discretisation`: Discretization parameters per mode

# Returns
- `TemporalHybridSymbolicModel`: The constructed symbolic model
"""
function Build_Timed_Hybrid_Automaton(hs::HybridSystem, growth_bounds::SVector{}, param_discretisation)
    # 1) Build symbolic models per mode
    abstract_systems = build_initial_symmodel_by_mode(hs, growth_bounds, param_discretisation)
    # 2) Build the global input map
    global_input_map = GlobalInputMap(abstract_systems, hs)
    # 3) Add intra-mode transitions
    transitions_list = []
    add_mode_transitions!(transitions_list, abstract_systems, global_input_map)
    # 4) Add switching transitions between modes
    add_switching_transitions!(transitions_list, hs, abstract_systems, global_input_map)
    # 5) Build the final automaton and state mappings
    int2aug_state, aug_state2int, autom = build_automaton(transitions_list)
    # 6) Extract symbolic and temporal models
    symmodels = [abs_sys[1] for abs_sys in abstract_systems]
    time_symbolic_models = [abs_sys[2] for abs_sys in abstract_systems]
    # 7) Return the complete temporal hybrid symbolic model
    return TemporalHybridSymbolicModel(
        symmodels,
        time_symbolic_models,
        int2aug_state,
        aug_state2int,
        autom,
        global_input_map
    )
end


# ================================================================
# Utility functions for transitions
# ================================================================

"""
    find_symbolic_state(symmodel, continuous_state)

Find the symbolic state index corresponding to a given continuous state.

# Arguments
- `symmodel`: The symbolic model
- `continuous_state`: The continuous state vector

# Returns
- `Int`: The symbolic state index (0 if not found)
"""
function find_symbolic_state(symmodel, continuous_state)
    state_idx = SY.get_abstract_state(symmodel, continuous_state)
    if isnothing(state_idx)
        return 0
    else
        return state_idx
    end
end

"""
    find_time_index(time_model, time_value)

Find the time index in the symbolic time model corresponding to a given time value.

# Arguments
- `time_model`: The symbolic time model
- `time_value`: The time value

# Returns
- `Int`: The time index (closest if not exact)
"""
function find_time_index(time_model, time_value)
    for (idx, tstep) in enumerate(time_model.tsteps)
        if abs(time_value - tstep) < 1e-10
            return idx
        end
    end
    # If no exact match, take the closest
    min_distance = Inf
    best_idx = 0
    for (idx, tstep) in enumerate(time_model.tsteps)
        distance = abs(time_value - tstep)
        if distance < min_distance
            min_distance = distance
            best_idx = idx
        end
    end
    return best_idx
end

# Pour le moment seul les guards basés sur des HyperRectangles sont supportés (à changer à terme).
"""
    extract_spatial_part(guard)

Extract the spatial part (all but last dimension) from a guard (assumed to be a HyperRectangle).

# Arguments
- `guard`: The guard set (should be a HyperRectangle)

# Returns
- `HyperRectangle`: The spatial part of the guard
"""
function extract_spatial_part(guard)
    if isa(guard, UT.HyperRectangle)
        return UT.HyperRectangle(guard.lb[1:end-1], guard.ub[1:end-1])
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
    if isa(guard, UT.HyperRectangle)
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
    indices = Int[]
    for (idx, tstep) in enumerate(time_model.tsteps)
        if t_min <= tstep <= t_max
            push!(indices, idx)
        end
    end
    return indices
end


end