# ------------------------------------------------
# Approximation kernels
# ------------------------------------------------

# One kernel for every over-approximation: `ST.input_cache` hoists whatever the
# approximation can precompute per input (for the uniform cell radius `r`), and
# `ST.reach_set` is the per-cell hot path. Approximation-specific formulas live in
# src/system/approximation/, not here.
function collect_abstract_transitions!(
    out::Vector{TransitionKey},
    metadata_out::Vector{Pair{TransitionKey, Any}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemOverApproximation;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    XMapping = get_state_mapping(symmodel)
    r = MP.get_h(MP.get_grid(XMapping)) / 2.0

    inputs = collect(enum_inputs(symmodel))
    input_data = Dict(
        u => begin
            concrete_input = get_concrete_input(symmodel, u)
            (concrete_input, ST.input_cache(concrete_system_approx, r, concrete_input))
        end for u in inputs
    )

    states = collect(enum_states(symmodel))
    concrete_elems = Dict(q => get_concrete_elem(symmodel, q) for q in states)

    workfun! = function (
        transbuf::Vector{TransitionKey},
        metadatabuf::Vector{Pair{TransitionKey, Any}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input, cache = input_data[abstract_input]
        concrete_elem = concrete_elems[abstract_state]
        reachable_set =
            ST.reach_set(concrete_system_approx, concrete_elem, concrete_input, cache)

        start_len = length(transbuf)

        allin = compute_abstract_transitions_from_set!(
            symmodel,
            reachable_set,
            abstract_state,
            abstract_input,
            transbuf,
        )

        allin || resize!(transbuf, start_len)
        return nothing
    end

    return _collect_transitions!(
        out,
        metadata_out,
        symmodel,
        workfun!;
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )
end

function collect_abstract_transitions!(
    out::Vector{TransitionKey},
    metadata_out::Vector{Pair{TransitionKey, Any}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeSystemUnderApproximation;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    under_approximation_map = ST.get_under_approximation_map(concrete_system_approx)

    inputs = collect(enum_inputs(symmodel))
    concrete_inputs = Dict(u => get_concrete_input(symmodel, u) for u in inputs)

    states = collect(enum_states(symmodel))
    concrete_states = Dict(q => get_concrete_state(symmodel, q) for q in states)
    concrete_elems = Dict(q => get_concrete_elem(symmodel, q) for q in states)

    workfun! = function (
        transbuf::Vector{TransitionKey},
        metadatabuf::Vector{Pair{TransitionKey, Any}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input = concrete_inputs[abstract_input]
        concrete_elem = concrete_elems[abstract_state]
        concrete_state = concrete_states[abstract_state]

        reachable_points = under_approximation_map(concrete_elem, concrete_input)

        start_len = length(transbuf)
        start_meta_len = length(metadatabuf)

        allin = compute_abstract_transitions_from_points!(
            symmodel,
            reachable_points,
            concrete_state,
            abstract_state,
            abstract_input,
            transbuf,
            metadatabuf,
        )

        if !allin
            resize!(transbuf, start_len)
            resize!(metadatabuf, start_meta_len)
        end

        return nothing
    end

    return _collect_transitions!(
        out,
        metadata_out,
        symmodel,
        workfun!;
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )
end

function collect_abstract_transitions!(
    out::Vector{TransitionKey},
    metadata_out::Vector{Pair{TransitionKey, Any}},
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx::ST.DiscreteTimeCenteredSimulation;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    system_map = ST.get_system_map(concrete_system_approx)

    inputs = collect(enum_inputs(symmodel))
    concrete_inputs = Dict(u => get_concrete_input(symmodel, u) for u in inputs)

    states = collect(enum_states(symmodel))
    concrete_states = Dict(q => get_concrete_state(symmodel, q) for q in states)

    workfun! = function (
        transbuf::Vector{TransitionKey},
        metadatabuf::Vector{Pair{TransitionKey, Any}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input = concrete_inputs[abstract_input]
        concrete_state = concrete_states[abstract_state]

        y = system_map(concrete_state, concrete_input)
        target = get_abstract_state(symmodel, y)

        if target !== nothing && is_allowed_state(symmodel, target)
            tr = (target, abstract_state, abstract_input)
            push!(transbuf, tr)

            if has_metadata(symmodel)
                push!(metadatabuf, tr => ConcreteTransitionSample(concrete_state, y))
            end
        end

        return nothing
    end

    return _collect_transitions!(
        out,
        metadata_out,
        symmodel,
        workfun!;
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )
end
