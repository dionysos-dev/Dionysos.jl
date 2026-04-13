# ------------------------------------------------
# Distributed Abstraction
# ------------------------------------------------

import Distributed

# ------------------------------------------------
# Worker-local persistent storage
# ------------------------------------------------

# These live once per Julia process. If Dionysos is loaded on each worker,
# each worker gets its own copy.
const _DIST_SYMMODEL = Ref{Any}(nothing)
const _DIST_APPROX = Ref{Any}(nothing)
const _DIST_READY = Ref(false)

# ------------------------------------------------
# Local partition model wrapper
# ------------------------------------------------

"""
    LocalGridBasedSymbolicModel

Wrapper around a global symbolic model that overrides only the source-domain
enumeration while keeping the same global state/input mappings and retained set.
"""
struct LocalGridBasedSymbolicModel{N, M, SM, XS} <: GridBasedSymbolicModel{N, M}
    parent::SM
    Xset_local::XS
end

LocalGridBasedSymbolicModel(
    symmodel::GridBasedSymbolicModel{N, M},
    Xset_local,
) where {N, M} = LocalGridBasedSymbolicModel{N, M, typeof(symmodel), typeof(Xset_local)}(
    symmodel,
    Xset_local,
)

get_state_mapping(sym::LocalGridBasedSymbolicModel) = get_state_mapping(sym.parent)
get_input_mapping(sym::LocalGridBasedSymbolicModel) = get_input_mapping(sym.parent)
get_source_domain(sym::LocalGridBasedSymbolicModel) = sym.Xset_local
get_retained_domain(sym::LocalGridBasedSymbolicModel) = get_retained_domain(sym.parent)
get_input_domain(sym::LocalGridBasedSymbolicModel) = get_input_domain(sym.parent)

# Forward everything else to parent
get_concrete_state(sym::LocalGridBasedSymbolicModel, q) = get_concrete_state(sym.parent, q)
get_concrete_elem(sym::LocalGridBasedSymbolicModel, q::Int) =
    get_concrete_elem(sym.parent, q)
get_concrete_input(sym::LocalGridBasedSymbolicModel, u) = get_concrete_input(sym.parent, u)
get_abstract_state(sym::LocalGridBasedSymbolicModel, x) = get_abstract_state(sym.parent, x)
add_transitions!(sym::LocalGridBasedSymbolicModel, trans) =
    add_transitions!(sym.parent, trans)

# ------------------------------------------------
# Partitioning
# ------------------------------------------------

"""
    partition_source_state_ids(symmodel, nparts; strategy=:roundrobin)

Partition the source abstract states into `nparts` lists of state ids.
This ships only integer ids to workers, not the whole symbolic model.
"""
function partition_source_state_ids(
    symmodel::GridBasedSymbolicModel,
    nparts::Int;
    strategy::Symbol = :roundrobin,
)
    nparts >= 1 || error("nparts must be >= 1")

    states = collect(enum_source_states(symmodel))
    isempty(states) && return [Int[] for _ in 1:nparts]

    parts = [Int[] for _ in 1:nparts]

    if strategy == :roundrobin
        for (k, q) in enumerate(states)
            push!(parts[mod1(k, nparts)], q)
        end
    elseif strategy == :contiguous
        n = length(states)
        for i in 1:nparts
            a = fld((i - 1) * n, nparts) + 1
            b = fld(i * n, nparts)
            if a <= b
                append!(parts[i], @view states[a:b])
            end
        end
    else
        error("Unknown partition strategy: $strategy")
    end

    return parts
end

"""
    assign_states_to_workers(parts, procs)

Assign the partition list `parts` to workers in `procs` using round-robin,
then merge the partitions assigned to the same worker into a single vector of
state ids. Returns a vector of pairs `(p, ids)`.
"""
function assign_states_to_workers(parts::Vector{Vector{Int}}, procs)
    isempty(procs) && return Tuple{Int, Vector{Int}}[]

    buckets = [Int[] for _ in 1:length(procs)]

    for (k, ids) in enumerate(parts)
        append!(buckets[mod1(k, length(procs))], ids)
    end

    return [(procs[i], buckets[i]) for i in eachindex(procs)]
end

# ------------------------------------------------
# Worker-local helpers
# ------------------------------------------------

"""
    _install_distributed_abstraction_data!(symmodel, concrete_system_approx)

Install the heavy read-only abstraction objects in worker-local storage.
This is meant to run on each worker process once.
"""
function _install_distributed_abstraction_data!(symmodel, concrete_system_approx)
    _DIST_SYMMODEL[] = symmodel
    _DIST_APPROX[] = concrete_system_approx
    _DIST_READY[] = true
    return nothing
end

"""
    _clear_distributed_abstraction_data!()

Clear worker-local persistent storage.
"""
function _clear_distributed_abstraction_data!()
    _DIST_SYMMODEL[] = nothing
    _DIST_APPROX[] = nothing
    _DIST_READY[] = false
    return nothing
end

"""
    _make_local_symmodel_from_ids(state_ids)

Construct a local wrapper on the current worker using the worker-local parent
symbolic model and the provided state ids.
"""
function _make_local_symmodel_from_ids(state_ids::Vector{Int})
    _DIST_READY[] || error("Distributed abstraction worker not initialized.")

    sym = _DIST_SYMMODEL[]
    Xmap = get_state_mapping(sym)
    Xset_local = MP.stateset_from_states(Xmap, state_ids)
    return LocalGridBasedSymbolicModel(sym, Xset_local)
end

"""
    _warmup_distributed_abstraction_worker!()

Cheap warm-up of the main abstraction code paths on the current worker.
This helps separate compilation/setup from actual benchmark timings.
"""
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

# ------------------------------------------------
# Result container
# ------------------------------------------------

struct DistributedAbstractionResult
    transitions::Vector{Tuple{Int, Int, Int}}
    n_source_states::Int
    n_transitions::Int
end

# ------------------------------------------------
# Worker kernel
# ------------------------------------------------

"""
    _run_local_partition_ids(state_ids; ...)

Run abstraction on the current worker for the source states listed in `state_ids`,
using the worker-local persistent symbolic model and concrete approximation.
"""
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

    transitions = collect_abstract_transitions(
        local_symmodel,
        local_approx;
        print_level = print_level,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
        state_filter = state_filter,
        state_input_filter = state_input_filter,
    )

    return DistributedAbstractionResult(transitions, length(state_ids), length(transitions))
end

# ------------------------------------------------
# Worker initialization API
# ------------------------------------------------

"""
    init_abstraction_workers!(procs, symmodel, concrete_system_approx; warmup=true)

Install heavy read-only data on each worker once. Optionally warm up the worker.
"""
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

"""
    clear_abstraction_workers!(procs)

Clear worker-local persistent storage.
"""
function clear_abstraction_workers!(procs = Distributed.workers())
    for p in procs
        Distributed.remotecall_wait(_clear_distributed_abstraction_data!, p)
    end
    return nothing
end

# ------------------------------------------------
# Public distributed mutating API
# ------------------------------------------------

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
        @info "Starting distributed abstraction" nprocs = length(procs) nparts = nparts partition_strategy =
            partition_strategy threaded_per_worker = threaded_per_worker master_nthreads =
            Threads.nthreads()
    end

    transitions = collect_abstract_transitions_distributed(
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

    # Install heavy read-only objects once per worker
    init_abstraction_workers!(
        procs,
        symmodel,
        concrete_system_approx;
        warmup = warmup_workers,
    )

    # First partition finely for balancing, then merge assigned partitions per worker
    parts = partition_source_state_ids(symmodel, nparts; strategy = partition_strategy)
    worker_assignments = assign_states_to_workers(parts, procs)

    futures = Vector{Distributed.Future}(undef, length(worker_assignments))

    for (i, (p, ids)) in enumerate(worker_assignments)
        futures[i] = Distributed.remotecall(
            _run_local_partition_ids,
            p,
            ids;
            print_level = (print_level <= 1 ? print_level : 1),
            update_interval = update_interval,
            progress_dt = progress_dt,
            threaded = threaded_per_worker,
            state_filter = state_filter,
            state_input_filter = state_input_filter,
        )
    end

    results = fetch.(futures)

    transitions = Tuple{Int, Int, Int}[]
    total_sources = 0
    total_transitions = 0

    for res in results
        append!(transitions, res.transitions)
        total_sources += res.n_source_states
        total_transitions += res.n_transitions
    end

    print_level >= 1 && @info(
        "Distributed abstraction finished",
        nparts = nparts,
        total_source_states = total_sources,
        total_transitions = total_transitions,
    )

    return transitions
end
