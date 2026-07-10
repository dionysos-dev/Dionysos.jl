# ------------------------------------------------------------
# Distributed Execution
# ------------------------------------------------------------

import Distributed

"""
    JuliaDistributedBackend(
        procs=nothing,
        nparts=nothing,
        partition_strategy=:roundrobin,
        threaded_per_worker=false,
    )

Distributed execution over Julia worker processes. Each worker receives its
share of the work in a single `remotecall` carrying the symbolic model and
the approximation explicitly — there is no per-worker global state to
install or clear, and workers JIT-compile in parallel inside their call.

# Parameters
- `procs`: worker IDs (defaults to `Distributed.workers()`).
- `nparts`: number of partitions (defaults to number of workers).
- `partition_strategy`: how to split states (`:roundrobin` or `:contiguous`).
- `threaded_per_worker`: enable threading inside each worker.
"""
struct JuliaDistributedBackend <: AbstractExecutionBackend
    procs::Union{Nothing, Vector{Int}}
    nparts::Union{Nothing, Int}
    partition_strategy::Symbol
    threaded_per_worker::Bool
end

JuliaDistributedBackend() = JuliaDistributedBackend(nothing, nothing, :roundrobin, false)

function compute_abstract_system_from_concrete_system!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    execution_backend::JuliaDistributedBackend;
    kwargs...,
)
    procs =
        execution_backend.procs === nothing ? Distributed.workers() :
        execution_backend.procs
    nparts =
        execution_backend.nparts === nothing ? max(length(procs), 1) :
        execution_backend.nparts

    return compute_abstract_system_distributed!(
        symmodel,
        concrete_system_approx;
        procs = procs,
        nparts = nparts,
        partition_strategy = execution_backend.partition_strategy,
        threaded_per_worker = execution_backend.threaded_per_worker,
        kwargs...,
    )
end

function assign_states_to_workers(parts::Vector{Vector{Int}}, procs)
    isempty(procs) && return Tuple{Int, Vector{Int}}[]

    buckets = [Int[] for _ in eachindex(procs)]

    for (k, ids) in enumerate(parts)
        append!(buckets[mod1(k, length(procs))], ids)
    end

    return [(procs[i], buckets[i]) for i in eachindex(procs)]
end

struct DistributedAbstractionResult
    transitions::Vector{TransitionKey}
    metadata_pairs::Vector{Pair{TransitionKey, Any}}
    n_source_states::Int
    n_transitions::Int
end

# Worker entry point: everything it needs arrives as arguments.
function _run_local_partition_ids(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    state_ids::Vector{Int};
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter = nothing,
    state_input_filter = nothing,
)
    local_symmodel = local_symmodel_from_state_ids(symmodel, state_ids)

    transitions, metadata_pairs = collect_abstract_transitions(
        local_symmodel,
        concrete_system_approx;
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )

    return DistributedAbstractionResult(
        transitions,
        metadata_pairs,
        length(state_ids),
        length(transitions),
    )
end

function compute_abstract_system_distributed!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    procs = Distributed.workers(),
    nparts::Int = max(length(procs), 1),
    partition_strategy::Symbol = :roundrobin,
    threaded_per_worker::Bool = false,
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    state_filter = nothing,
    state_input_filter = nothing,
)
    if print_level >= 1
        @info(
            "Starting distributed abstraction",
            nprocs = length(procs),
            nparts = nparts,
            partition_strategy = partition_strategy,
            threaded_per_worker = threaded_per_worker,
            master_nthreads = Threads.nthreads(),
        )
    end

    transitions, metadata_pairs = collect_abstract_transitions_distributed(
        symmodel,
        concrete_system_approx;
        procs = procs,
        nparts = nparts,
        partition_strategy = partition_strategy,
        threaded_per_worker = threaded_per_worker,
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )

    isempty(transitions) || add_transitions!(symmodel, transitions)
    add_metadata_pairs!(symmodel, metadata_pairs)
    return symmodel
end

function collect_abstract_transitions_distributed(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    procs = Distributed.workers(),
    nparts::Int = max(length(procs), 1),
    partition_strategy::Symbol = :roundrobin,
    threaded_per_worker::Bool = false,
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    state_filter = nothing,
    state_input_filter = nothing,
)
    nparts >= 1 || error("nparts must be >= 1")

    if isempty(procs)
        return collect_abstract_transitions(
            symmodel,
            concrete_system_approx;
            print_level = print_level,
            update_interval = update_interval,
            progress_dt = progress_dt,
            threaded = threaded_per_worker,
            state_filter = state_filter,
            state_input_filter = state_input_filter,
        )
    end

    parts = partition_source_state_ids(symmodel, nparts; strategy = partition_strategy)
    worker_assignments = assign_states_to_workers(parts, procs)

    if print_level >= 1
        ninputs = length(collect(enum_inputs(symmodel)))
        max_work = maximum(length(ids) * ninputs for (_, ids) in worker_assignments)
        @info "Max worker workload" max_state_input_checks = max_work
    end

    futures = Vector{Distributed.Future}(undef, length(worker_assignments))

    for (i, (p, ids)) in enumerate(worker_assignments)
        futures[i] = Distributed.remotecall(
            _run_local_partition_ids,
            p,
            symmodel,
            concrete_system_approx,
            ids;
            print_level = min(print_level, 1),
            update_interval = update_interval,
            progress_dt = progress_dt,
            threaded = threaded_per_worker,
            state_filter = state_filter,
            state_input_filter = state_input_filter,
        )
    end

    results = fetch.(futures)

    transitions = TransitionKey[]
    metadata_pairs = Pair{TransitionKey, Any}[]

    total_sources = 0
    total_transitions = 0

    for res in results
        append!(transitions, res.transitions)
        append!(metadata_pairs, res.metadata_pairs)

        total_sources += res.n_source_states
        total_transitions += res.n_transitions
    end

    print_level >= 1 && @info(
        "Distributed abstraction finished",
        nparts = nparts,
        total_source_states = total_sources,
        total_transitions = total_transitions,
    )

    return transitions, metadata_pairs
end
