# ------------------------------------------------
# Abstraction orchestration (backend-agnostic)
# ------------------------------------------------

function compute_abstract_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemApproximation;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    trans, metadata_pairs = collect_abstract_transitions(
        abstract_system,
        concrete_system_approx;
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )

    isempty(trans) || add_transitions!(abstract_system, trans)

    # Deduplicate / compress automaton after all transitions are inserted.
    finalize!(get_automaton(abstract_system))

    add_metadata_pairs!(abstract_system, metadata_pairs)

    return abstract_system
end

function collect_abstract_transitions(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    trans = TransitionKey[]
    metadata_pairs = Pair{TransitionKey, Any}[]

    collect_abstract_transitions!(
        trans,
        metadata_pairs,
        symmodel,
        concrete_system_approx;
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )

    return trans, metadata_pairs
end

# ------------------------------------------------
# Low-level transition helpers
# ------------------------------------------------

function compute_abstract_transitions_from_set!(
    symmodel::GridBasedSymbolicModel,
    reachable_set::LazySets.LazySet,
    abstract_state::Int,
    abstract_input::Int,
    translist::Vector{TransitionKey},
)
    targets, allin = get_states_from_set_strict(symmodel, reachable_set, MP.OUTER)
    allin || return false

    for target in targets
        is_allowed_state(symmodel, target) || return false
    end

    for target in targets
        push!(translist, (target, abstract_state, abstract_input))
    end

    return true
end

function compute_abstract_transitions_from_points!(
    symmodel::GridBasedSymbolicModel,
    reachable_points,
    concrete_source,
    abstract_state::Int,
    abstract_input::Int,
    translist::Vector{TransitionKey},
    metadatalist::Vector{Pair{TransitionKey, Any}},
)
    start_len = length(translist)
    start_meta_len = length(metadatalist)

    seen_targets = Set{Int}()

    for y in reachable_points
        target = get_abstract_state(symmodel, y)

        if target === nothing || !is_allowed_state(symmodel, target)
            resize!(translist, start_len)
            resize!(metadatalist, start_meta_len)
            return false
        end

        target in seen_targets && continue
        push!(seen_targets, target)

        tr = (target, abstract_state, abstract_input)
        push!(translist, tr)

        if has_metadata(symmodel)
            push!(metadatalist, tr => ConcreteTransitionSample(concrete_source, y))
        end
    end

    return true
end
