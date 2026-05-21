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
        warmup_workers=true,
    )

Distributed execution over Julia worker processes.

# Parameters
- `procs`: worker IDs (defaults to `Distributed.workers()`).
- `nparts`: number of partitions (defaults to number of workers).
- `partition_strategy`: how to split states (`:roundrobin` or `:contiguous`).
- `threaded_per_worker`: enable threading inside each worker.
- `warmup_workers`: run a small warm-up to reduce compilation overhead.
"""
struct JuliaDistributedBackend <: AbstractExecutionBackend
    procs::Union{Nothing, Vector{Int}}
    nparts::Union{Nothing, Int}
    partition_strategy::Symbol
    threaded_per_worker::Bool
    warmup_workers::Bool
end

JuliaDistributedBackend() =
    JuliaDistributedBackend(nothing, nothing, :roundrobin, false, true)

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
        warmup_workers = execution_backend.warmup_workers,
        kwargs...,
    )
end

const _DIST_SYMMODEL = Ref{Any}(nothing)
const _DIST_APPROX = Ref{Any}(nothing)
const _DIST_READY = Ref(false)

function assign_states_to_workers(parts::Vector{Vector{Int}}, procs)
    isempty(procs) && return Tuple{Int, Vector{Int}}[]

    buckets = [Int[] for _ in eachindex(procs)]

    for (k, ids) in enumerate(parts)
        append!(buckets[mod1(k, length(procs))], ids)
    end

    return [(procs[i], buckets[i]) for i in eachindex(procs)]
end

function _install_distributed_abstraction_data!(symmodel, concrete_system_approx)
    _DIST_SYMMODEL[] = symmodel
    _DIST_APPROX[] = concrete_system_approx
    _DIST_READY[] = true
    return nothing
end

function _clear_distributed_abstraction_data!()
    _DIST_SYMMODEL[] = nothing
    _DIST_APPROX[] = nothing
    _DIST_READY[] = false
    return nothing
end

function _make_local_symmodel_from_ids(state_ids::Vector{Int})
    _DIST_READY[] || error("Distributed abstraction worker not initialized.")

    sym = _DIST_SYMMODEL[]
    return local_symmodel_from_state_ids(sym, state_ids)
end

function _warmup_distributed_abstraction_worker!()
    _DIST_READY[] || error("Distributed abstraction worker not initialized.")

    sym = _DIST_SYMMODEL[]
    approx = _DIST_APPROX[]

    states = collect(enum_source_states(sym))
    inputs = collect(enum_inputs(sym))

    isempty(states) && return nothing
    isempty(inputs) && return nothing

    q = first(states)
    u = first(inputs)

    x = get_concrete_state(sym, q)
    cu = get_concrete_input(sym, u)

    if approx isa ST.DiscreteTimeCenteredSimulation
        f = ST.get_system_map(approx)
        f(x, cu)

    elseif approx isa ST.DiscreteTimeGrowthBound
        f = ST.get_system_map(approx)
        g = approx.growthbound_map
        f(x, cu)

        XMapping = get_state_mapping(sym)
        r = MP.get_h(MP.get_grid(XMapping)) / 2.0
        g(r, cu)

    elseif approx isa ST.DiscreteTimeSystemOverApproximation
        f = ST.get_over_approximation_map(approx)
        elem = get_concrete_elem(sym, q)
        f(elem, cu)

    elseif approx isa ST.DiscreteTimeSystemUnderApproximation
        f = ST.get_under_approximation_map(approx)
        elem = get_concrete_elem(sym, q)
        f(elem, cu)

    elseif approx isa ST.DiscreteTimeLinearized
        XMapping = get_state_mapping(sym)
        N = MP.get_dim(XMapping)
        r = MP.get_h(MP.get_grid(XMapping)) / 2.0
        H = SMatrix{N, N}(LA.I) .* r
        lmap = approx.linsys_map
        emap = approx.error_map
        lmap(x, H, cu)
        e = LA.norm(r, Inf)
        emap(e, cu)
    end

    return nothing
end

struct DistributedAbstractionResult
    transitions::Vector{TransitionKey}
    metadata_pairs::Vector{Pair{TransitionKey, Any}}
    n_source_states::Int
    n_transitions::Int
end

function _run_local_partition_ids(
    state_ids::Vector{Int};
    print_level::Int = 0,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
    state_filter = nothing,
    state_input_filter = nothing,
)
    local_symmodel = _make_local_symmodel_from_ids(state_ids)
    local_approx = _DIST_APPROX[]

    transitions, metadata_pairs = collect_abstract_transitions(
        local_symmodel,
        local_approx;
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

function init_abstraction_workers!(
    procs,
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    warmup::Bool = true,
)
    for p in procs
        Distributed.remotecall_wait(
            _install_distributed_abstraction_data!,
            p,
            symmodel,
            concrete_system_approx,
        )
    end

    if warmup
        for p in procs
            Distributed.remotecall_wait(_warmup_distributed_abstraction_worker!, p)
        end
    end

    return nothing
end

function clear_abstraction_workers!(procs = Distributed.workers())
    for p in procs
        Distributed.remotecall_wait(_clear_distributed_abstraction_data!, p)
    end

    return nothing
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
    warmup_workers::Bool = true,
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
        warmup_workers = warmup_workers,
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
    warmup_workers::Bool = true,
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

    init_abstraction_workers!(
        procs,
        symmodel,
        concrete_system_approx;
        warmup = warmup_workers,
    )

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
