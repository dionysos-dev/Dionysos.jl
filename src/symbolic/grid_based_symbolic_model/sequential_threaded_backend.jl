# ------------------------------------------------
# Sequential and Threaded Execution
# ------------------------------------------------

"""
    SequentialBackend()

Sequential execution (no parallelism).
"""
struct SequentialBackend <: AbstractExecutionBackend end

"""
    ThreadedBackend(progress_dt=0.2)

Multithreaded execution using all available Julia threads.
"""
struct ThreadedBackend <: AbstractExecutionBackend
    progress_dt::Float64
end

ThreadedBackend() = ThreadedBackend(0.2)

function compute_abstract_system_from_concrete_system!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    ::SequentialBackend;
    kwargs...,
)
    compute_abstract_system!(symmodel, concrete_system_approx; threaded = false, kwargs...)
    return symmodel
end

function compute_abstract_system_from_concrete_system!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    execution_backend::ThreadedBackend;
    kwargs...,
)
    compute_abstract_system!(
        symmodel,
        concrete_system_approx;
        threaded = true,
        progress_dt = execution_backend.progress_dt,
        kwargs...,
    )
    return symmodel
end

function compute_abstract_system!(
    abstract_system::GridBasedSymbolicModel,
    concrete_system_approx::Union{
        ST.DiscreteTimeSystemOverApproximation,
        ST.DiscreteTimeGrowthBound,
        ST.DiscreteTimeLinearized,
        ST.DiscreteTimeSystemUnderApproximation,
        ST.DiscreteTimeCenteredSimulation,
    };
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
    ST.finalize!(get_automaton(abstract_system))

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
# Generic helpers
# ------------------------------------------------

function _collect_transitions!(
    out::Vector{TransitionKey},
    metadata_out::Vector{Pair{TransitionKey, Any}},
    symmodel::GridBasedSymbolicModel,
    workfun!;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    if !threaded || Threads.nthreads() == 1
        return _collect_transitions_sequential!(
            out,
            metadata_out,
            symmodel,
            workfun!;
            print_level = print_level,
            update_interval = update_interval,
            state_filter = state_filter,
            state_input_filter = state_input_filter,
        )
    else
        return _collect_transitions_threaded!(
            out,
            metadata_out,
            symmodel,
            workfun!;
            print_level = print_level,
            update_interval = update_interval,
            progress_dt = progress_dt,
            state_filter = state_filter,
            state_input_filter = state_input_filter,
        )
    end
end

function _collect_transitions_sequential!(
    out::Vector{TransitionKey},
    metadata_out::Vector{Pair{TransitionKey, Any}},
    symmodel::GridBasedSymbolicModel,
    workfun!;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    print_level >= 1 && @info("Starting sequential abstraction")

    inputs = collect(enum_inputs(symmodel))
    states = filtered_source_states(symmodel, state_filter)

    total_work = length(inputs) * length(states)
    total_updates = max(div(total_work, max(1, update_interval)), 1)
    progress = print_level == 2 ? ProgressMeter.Progress(total_updates) : nothing

    localbuf = TransitionKey[]
    localmeta = Pair{TransitionKey, Any}[]
    count = 0

    for abstract_input in inputs
        for abstract_state in states
            empty!(localbuf)
            empty!(localmeta)

            if state_input_filter === nothing || _keep_state_input(
                symmodel,
                abstract_state,
                abstract_input,
                state_input_filter,
            )
                workfun!(localbuf, localmeta, abstract_state, abstract_input)
                append!(out, localbuf)
                append!(metadata_out, localmeta)
            end

            count += 1

            if print_level >= 1 && count % update_interval == 0
                if print_level == 1
                    @info "Sequential abstraction progress" done = count total = total_work
                elseif print_level == 2
                    ProgressMeter.next!(progress)
                end
            end
        end
    end

    print_level == 2 && ProgressMeter.finish!(progress)

    print_level >= 1 && @info(
        "Finished sequential abstraction",
        total_work = total_work,
        ntransitions = length(out),
    )

    return out
end

function _collect_transitions_threaded!(
    out::Vector{TransitionKey},
    metadata_out::Vector{Pair{TransitionKey, Any}},
    symmodel::GridBasedSymbolicModel,
    workfun!;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    print_level >= 1 &&
        @info("Starting threaded abstraction", nthreads = Threads.nthreads())

    inputs = collect(enum_inputs(symmodel))
    states = filtered_source_states(symmodel, state_filter)

    total_work = length(inputs) * length(states)

    nthreads = Threads.nthreads()
    nbuffers = Threads.maxthreadid()

    transitions_by_thread = [TransitionKey[] for _ in 1:nbuffers]
    metadata_by_thread = [Pair{TransitionKey, Any}[] for _ in 1:nbuffers]
    local_done = fill(0, nbuffers)

    progress_dt_ns = Int(round(progress_dt * 1e9))
    prog = print_level == 2 ? ProgressMeter.Progress(total_work) : nothing
    global_done = Threads.Atomic{Int}(0)
    last_t = time_ns()

    Threads.@threads for linear_idx in 1:total_work
        tid = Threads.threadid()

        local_transitions = transitions_by_thread[tid]
        local_metadata = metadata_by_thread[tid]

        input_idx = ((linear_idx - 1) ÷ length(states)) + 1
        state_idx = ((linear_idx - 1) % length(states)) + 1

        abstract_input = @inbounds inputs[input_idx]
        abstract_state = @inbounds states[state_idx]

        if state_input_filter === nothing ||
           _keep_state_input(symmodel, abstract_state, abstract_input, state_input_filter)
            workfun!(local_transitions, local_metadata, abstract_state, abstract_input)
        end

        local_done[tid] += 1

        if local_done[tid] >= update_interval
            Threads.atomic_add!(global_done, local_done[tid])
            local_done[tid] = 0
        end

        if print_level >= 1 && tid == 1 && print_level == 2
            t = time_ns()
            if t - last_t >= progress_dt_ns
                ProgressMeter.update!(prog, global_done[] + local_done[1])
                last_t = t
            end
        end
    end

    for tid in 1:nbuffers
        if local_done[tid] > 0
            Threads.atomic_add!(global_done, local_done[tid])
        end
    end

    for tid in 1:nbuffers
        isempty(transitions_by_thread[tid]) || append!(out, transitions_by_thread[tid])
        isempty(metadata_by_thread[tid]) || append!(metadata_out, metadata_by_thread[tid])
    end

    if print_level == 2
        ProgressMeter.update!(prog, global_done[])
        ProgressMeter.finish!(prog)
    end

    print_level >= 1 && @info(
        "Finished threaded abstraction",
        total_work = total_work,
        nthreads = nthreads,
        ntransitions = length(out),
    )

    return out
end

# ------------------------------------------------
# Low-level transition helpers
# ------------------------------------------------

function compute_abstract_transitions_from_rectangle!(
    symmodel::GridBasedSymbolicModel,
    reachable_set::UT.HyperRectangle,
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

# ------------------------------------------------
# Approximation kernels
# ------------------------------------------------

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
    compute_reachable_set = ST.get_over_approximation_map(concrete_system_approx)

    inputs = collect(enum_inputs(symmodel))
    concrete_inputs = Dict(u => get_concrete_input(symmodel, u) for u in inputs)

    states = collect(enum_states(symmodel))
    concrete_elems = Dict(q => get_concrete_elem(symmodel, q) for q in states)

    workfun! = function (
        transbuf::Vector{TransitionKey},
        metadatabuf::Vector{Pair{TransitionKey, Any}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input = concrete_inputs[abstract_input]
        concrete_elem = concrete_elems[abstract_state]
        reachable_set = compute_reachable_set(concrete_elem, concrete_input)

        start_len = length(transbuf)

        allin = compute_abstract_transitions_from_rectangle!(
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
    concrete_system_approx::ST.DiscreteTimeGrowthBound;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    growthbound_map = concrete_system_approx.growthbound_map
    system_map = ST.get_system_map(concrete_system_approx)

    XMapping = get_state_mapping(symmodel)
    r = MP.get_h(MP.get_grid(XMapping)) / 2.0

    inputs = collect(enum_inputs(symmodel))
    input_data = Dict{Int, Tuple{Any, Any}}()

    for abstract_input in inputs
        concrete_input = get_concrete_input(symmodel, abstract_input)
        Fr = growthbound_map(r, concrete_input)
        input_data[abstract_input] = (concrete_input, Fr)
    end

    states = collect(enum_states(symmodel))
    concrete_states = Dict(q => get_concrete_state(symmodel, q) for q in states)

    workfun! = function (
        transbuf::Vector{TransitionKey},
        metadatabuf::Vector{Pair{TransitionKey, Any}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input, Fr = input_data[abstract_input]
        concrete_state = concrete_states[abstract_state]

        Fx = system_map(concrete_state, concrete_input)
        reachable_set = UT.HyperRectangle(Fx - Fr, Fx + Fr)

        start_len = length(transbuf)

        allin = compute_abstract_transitions_from_rectangle!(
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
    concrete_system_approx::ST.DiscreteTimeLinearized;
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter::Union{Nothing, Function} = nothing,
    state_input_filter::Union{Nothing, Function} = nothing,
)
    XMapping = get_state_mapping(symmodel)
    N = MP.get_dim(XMapping)
    r = MP.get_h(MP.get_grid(XMapping)) / 2.0
    _H_ = SMatrix{N, N}(LA.I) .* r
    _ONE_ = ones(SVector{N})
    e = LA.norm(r, Inf)

    error_map = concrete_system_approx.error_map
    linsys_map = concrete_system_approx.linsys_map

    inputs = collect(enum_inputs(symmodel))
    input_data = Dict{Int, Tuple{Any, Any, Any}}()

    for abstract_input in inputs
        concrete_input = get_concrete_input(symmodel, abstract_input)
        Fe = error_map(e, concrete_input)
        Fr = r .+ Fe
        input_data[abstract_input] = (concrete_input, Fe, Fr)
    end

    states = collect(enum_states(symmodel))
    concrete_states = Dict(q => get_concrete_state(symmodel, q) for q in states)

    workfun! = function (
        transbuf::Vector{TransitionKey},
        metadatabuf::Vector{Pair{TransitionKey, Any}},
        abstract_state::Int,
        abstract_input::Int,
    )
        concrete_input, Fe, Fr = input_data[abstract_input]
        concrete_state = concrete_states[abstract_state]

        Fx, DFx = linsys_map(concrete_state, _H_, concrete_input)

        A = inv(DFx)
        b = abs.(A) * Fr .+ 1.0

        rad = abs.(DFx) * _ONE_ .+ Fe
        reachable_set = UT.HyperRectangle(Fx - rad, Fx + rad)

        start_len = length(transbuf)

        allin = compute_abstract_transitions_from_rectangle!(
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
