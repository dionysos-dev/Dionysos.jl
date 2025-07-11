abstract type AbstractAutomatonList{N, M} <: HybridSystems.AbstractAutomaton end

# === Required Interface ===
function get_n_state(autom::AbstractAutomatonList{N, M}) where {N, M} end
function get_n_input(autom::AbstractAutomatonList{N, M}) where {N, M} end
function enum_transitions(autom::AbstractAutomatonList{N, M}) where {N, M} end
function add_transition!(
    autom::AbstractAutomatonList{N, M},
    source::Int,
    target::Int,
    symbol::Int,
) where {N, M} end
function pre(autom::AbstractAutomatonList{N, M}, target::Int) where {N, M} end
function post(autom::AbstractAutomatonList{N, M}, source::Int, symbol::Int) where {N, M} end
function Base.empty!(autom::AbstractAutomatonList{N, M}) where {N, M} end
function add_state!(autom::AbstractAutomatonList{N, M}) where {N, M} end

# === Common Default Implementations ===
enum_states(autom::AbstractAutomatonList) = 1:get_n_state(autom)
enum_inputs(autom::AbstractAutomatonList) = 1:get_n_input(autom)

function HybridSystems.ntransitions(autom::AbstractAutomatonList{N, M}) where {N, M}
    return length(enum_transitions(autom))
end

function add_transitions!(autom::AbstractAutomatonList{N, M}, translist) where {N, M}
    for (q′, q, u) in translist
        add_transition!(autom, q, q′, u)
    end
end

function is_deterministic(autom::AbstractAutomatonList{N, M}) where {N, M}
    seen = Dict{Tuple{Int, Int}, Int}()
    for (q′, q, u) in enum_transitions(autom)
        key = (q, u)
        seen[key] = get(seen, key, 0) + 1
        if seen[key] > 1
            return false
        end
    end
    return true
end

###################################################
################# Optimal control #################
###################################################

# Computes a controller to reach a target set in a non-deterministic automaton with a worst-case cost model.
# Assumes that states and inputs of the automaton are [1, 2, ..., nstates] and [1, 2, ..., nsymbols]
function compute_worst_case_cost_controller(
    autom::AbstractAutomatonList,
    target_set;
    initial_set = enum_states(autom),
    sparse_input = false,
    cost_function = nothing,
    ControllerConstructor::Function = () -> ST.SymbolicControllerList(),
)
    if cost_function === nothing
        abstract_controller, controllable_set, uncontrollable_set, value_fun_tab =
            compute_worst_case_uniform_cost_controller(
                autom,
                target_set;
                initial_set = initial_set,
                sparse_input = sparse_input,
                ControllerConstructor = ControllerConstructor,
            )
    else
        if is_deterministic(autom)
            abstract_controller, controllable_set, uncontrollable_set, value_fun_tab =
                compute_dijkstra_controller(
                    autom,
                    target_set;
                    initial_set = initial_set,
                    cost_function = cost_function,
                    ControllerConstructor = ControllerConstructor,
                )
        else
            abstract_controller, controllable_set, uncontrollable_set, value_fun_tab =
                compute_worst_case_non_uniform_cost_controller(
                    autom,
                    target_set;
                    initial_set = initial_set,
                    sparse_input = sparse_input,
                    cost_function = cost_function,
                    ControllerConstructor = ControllerConstructor,
                )
        end
    end

    return abstract_controller, controllable_set, uncontrollable_set, value_fun_tab
end
function compute_worst_case_uniform_cost_controller(
    autom::AbstractAutomatonList,
    target_set;
    initial_set = enum_states(autom),
    sparse_input = false,
    ControllerConstructor::Function = () -> ST.SymbolicControllerList(),
)
    abstract_controller = ControllerConstructor()
    stateset,
    initset,
    controllable_set,
    num_targets_unreachable,
    current_targets,
    next_targets,
    value_fun_tab = _data(autom, initial_set, target_set, sparse_input)

    success, value_fun_tab = _compute_controller_reach!(
        abstract_controller,
        autom,
        initset,
        controllable_set,
        num_targets_unreachable,
        current_targets,
        next_targets,
        value_fun_tab,
    )

    uncontrollable_set = setdiff(stateset, controllable_set)

    return abstract_controller, controllable_set, uncontrollable_set, value_fun_tab
end

function increase_counter!(counter::Array{Int, 2}, source::Int, symbol::Int)
    return counter[source, symbol] += 1
end
function increase_counter!(counter::Dict{Tuple{Int, Int}, Int}, source::Int, symbol::Int)
    key = (source, symbol)
    return counter[key] = get(counter, key, 0) + 1
end

function decrease_counter!(counter::Array{Int, 2}, source::Int, symbol::Int)
    counter[source, symbol] -= 1
    return counter[source, symbol]
end
function decrease_counter!(counter::Dict{Tuple{Int, Int}, Int}, source::Int, symbol::Int)
    key = (source, symbol)
    counter[key] = get(counter, key, 0) - 1
    return counter[key]
end

function _compute_num_targets_unreachable(counter, autom)
    for target in enum_states(autom)
        for (source, symbol) in pre(autom, target)
            increase_counter!(counter, source, symbol)
        end
    end
end

function _data(autom, initlist, targetlist, sparse_input::Bool)
    if sparse_input
        num_targets_unreachable = Dict{Tuple{Int, Int}, Int}()
    else
        num_targets_unreachable = zeros(Int, get_n_state(autom), get_n_input(autom))
    end

    _compute_num_targets_unreachable(num_targets_unreachable, autom)

    stateset = BitSet(enum_states(autom))
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    value_fun_tab = fill(Inf, get_n_state(autom)) # Inf = uncontrollable by default

    return stateset,
    initset,
    targetset,
    num_targets_unreachable,
    current_targets,
    next_targets,
    value_fun_tab
end

function _compute_controller_reach!(
    contr,
    autom,
    init_set,
    target_set,
    counter,
    current_targets,
    next_targets,
    value_fun_tab,
)::Tuple{Bool, Vector{Float64}}
    num_init_unreachable = length(init_set)

    step = 0
    for s in current_targets
        value_fun_tab[s] = step
    end

    while !isempty(current_targets) && !iszero(num_init_unreachable)
        empty!(next_targets)
        step += 1

        for target in current_targets
            for (source, symbol) in pre(autom, target)
                if !(source in target_set) &&
                   iszero(decrease_counter!(counter, source, symbol))
                    push!(target_set, source)
                    push!(next_targets, source)
                    ST.add_control!(contr, source, symbol)
                    value_fun_tab[source] = step

                    if source in init_set
                        num_init_unreachable -= 1
                    end
                end
            end
        end

        current_targets, next_targets = next_targets, current_targets
    end

    return iszero(num_init_unreachable), value_fun_tab
end

function compute_worst_case_non_uniform_cost_controller(
    autom::AbstractAutomatonList,
    target_set;
    initial_set = enum_states(autom),
    sparse_input = false,
    cost_function = (q, u) -> 1.0,
    ControllerConstructor::Function = () -> ST.SymbolicControllerList(),
)
    abstract_controller = ControllerConstructor()
    stateset, initset, targetset, counter, current_targets, next_targets, value_fun_tab =
        _data(autom, initial_set, target_set, sparse_input)

    num_init_unreachable = length(initset)

    # Set initial costs to 0
    for s in current_targets
        value_fun_tab[s] = 0.0
    end

    while !isempty(current_targets) && num_init_unreachable > 0
        empty!(next_targets)

        for target in current_targets
            for (source, symbol) in pre(autom, target)
                if source in targetset
                    continue
                end

                if decrease_counter!(counter, source, symbol) == 0
                    # All successors of (source, symbol) are in target set now
                    # Compute worst-case cost
                    successors = post(autom, source, symbol)

                    if isempty(successors)
                        continue
                    end
                    worst = maximum(value_fun_tab[q′] for q′ in successors)
                    total_cost = cost_function(source, symbol) + worst

                    if total_cost < value_fun_tab[source]
                        value_fun_tab[source] = total_cost
                        ST.set_control!(abstract_controller, source, symbol)
                        push!(targetset, source)
                        push!(next_targets, source)

                        if source in initset
                            num_init_unreachable -= 1
                        end
                    end
                end
            end
        end

        current_targets, next_targets = next_targets, current_targets
    end

    uncontrollable_set = setdiff(stateset, targetset)
    return abstract_controller, targetset, uncontrollable_set, value_fun_tab
end

using DataStructures  # For PriorityQueue

function compute_dijkstra_controller(
    autom::AbstractAutomatonList,
    target_set;
    initial_set = enum_states(autom),
    cost_function = (q, u) -> 1.0,
    ControllerConstructor::Function = () -> ST.SymbolicControllerList(),
)
    abstract_controller = ControllerConstructor()
    stateset = enum_states(autom)

    # Initialize value function: ∞ for all, 0 for target states
    value_fun_tab = Dict(q => Inf for q in stateset)
    for q in target_set
        value_fun_tab[q] = 0.0
    end

    # Priority queue keyed by value (cost)
    pq = PriorityQueue{Int, Float64}()
    for q in target_set
        enqueue!(pq, q, 0.0)
    end

    while !isempty(pq)
        (target, cost_to_target) = dequeue_pair!(pq)
        for (source, symbol) in pre(autom, target)
            new_cost = cost_function(source, symbol) + cost_to_target

            if new_cost < value_fun_tab[source]
                value_fun_tab[source] = new_cost
                enqueue!(pq, source, new_cost)
                ST.set_control!(abstract_controller, source, symbol)
            end
        end
    end

    reachable = Set(k for (k, v) in value_fun_tab if isfinite(v))
    uncontrollable = setdiff(initial_set, reachable)

    return abstract_controller, reachable, uncontrollable, value_fun_tab
end

###################################################
################# Safety control ##################
###################################################

function compute_largest_invariant_set(
    autom::AbstractAutomatonList,
    safelist;
    ControllerConstructor::Function = () -> ST.SymbolicControllerList(),
)
    contr = ControllerConstructor()
    nstates = get_n_state(autom)
    nsymbols = get_n_input(autom)
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]

    _compute_pairstable(pairstable, autom)
    nsymbolslist = sum(pairstable; dims = 2)

    # Remove unsafe states
    safeset = Set(safelist)
    for source in safeset
        if nsymbolslist[source] == 0
            delete!(safeset, source)
        end
    end

    unsafeset = Set(1:nstates)
    setdiff!(unsafeset, safeset)

    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end
    nextunsafeset = Set{Int}()

    # Iterate until convergence
    while true
        for target in unsafeset
            for soursymb in pre(autom, target)
                if pairstable[soursymb[1], soursymb[2]]
                    pairstable[soursymb[1], soursymb[2]] = false
                    nsymbolslist[soursymb[1]] -= 1
                    if nsymbolslist[soursymb[1]] == 0
                        push!(nextunsafeset, soursymb[1])
                    end
                end
            end
        end

        if isempty(nextunsafeset)
            break
        end

        setdiff!(safeset, nextunsafeset)
        unsafeset, nextunsafeset = nextunsafeset, unsafeset
        empty!(nextunsafeset)
    end

    # Populate controller
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                ST.add_control!(contr, source, symbol)
            end
        end
    end
    unsafeset = setdiff(Set(safelist), safeset)
    return contr, safeset, unsafeset
end

function _compute_pairstable(pairstable, autom)
    for target in enum_states(autom)
        for soursymb in pre(autom, target)
            pairstable[soursymb[1], soursymb[2]] = true
        end
    end
end
