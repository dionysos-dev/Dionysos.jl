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

# ------------------------------------------------
# Local execution primitive (sequential / threaded)
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
