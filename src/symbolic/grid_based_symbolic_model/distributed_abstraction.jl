# ------------------------------------------------
# Disibuted Abstraction
# ------------------------------------------------
import Distributed

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

# forward everything else to parent
get_concrete_state(sym::LocalGridBasedSymbolicModel, q) = get_concrete_state(sym.parent, q)
get_concrete_elem(sym::LocalGridBasedSymbolicModel, q) = get_concrete_elem(sym.parent, q)
get_concrete_input(sym::LocalGridBasedSymbolicModel, u) = get_concrete_input(sym.parent, u)
get_abstract_state(sym::LocalGridBasedSymbolicModel, x) = get_abstract_state(sym.parent, x)
add_transitions!(sym::LocalGridBasedSymbolicModel, trans) =
    add_transitions!(sym.parent, trans)

# ------------------------------------------------------------------
# Partitioning
# ------------------------------------------------------------------

function partition_source_states(
    symmodel::GridBasedSymbolicModel,
    nparts::Int;
    strategy::Symbol = :roundrobin,
)
    nparts >= 1 || error("nparts must be >= 1")

    states = collect(enum_source_states(symmodel))
    isempty(states) && return [ExplicitStateSubset(Int[]) for _ in 1:nparts]

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

    Xmap = get_state_mapping(symmodel)
    return [MP.stateset_from_states(Xmap, p) for p in parts]
end

# ------------------------------------------------------------------
# Local worker kernel
# ------------------------------------------------------------------

struct DistributedAbstractionResult
    transitions::Vector{Tuple{Int, Int, Int}}
    n_source_states::Int
    n_transitions::Int
end

function _run_local_partition(
    local_symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    verbose::Bool = false,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
    threaded::Bool = false,
)
    transitions = collect_abstract_transitions(
        local_symmodel,
        concrete_system_approx;
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
        threaded = threaded,
    )

    n_source_states = length(collect(enum_source_states(local_symmodel)))

    return DistributedAbstractionResult(transitions, n_source_states, length(transitions))
end

# ------------------------------------------------------------------
# Public distributed mutating API
# ------------------------------------------------------------------

function compute_abstract_system_distributed!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    procs = Distributed.workers(),
    nparts::Int = length(procs),
    partition_strategy::Symbol = :roundrobin,
    threaded_per_worker::Bool = false,
    verbose::Bool = true,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
)
    if verbose
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
        verbose = verbose,
        update_interval = update_interval,
        progress_dt = progress_dt,
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
    verbose::Bool = true,
    update_interval::Int = Int(1e5),
    progress_dt::Float64 = 0.2,
)
    nparts >= 1 || error("nparts must be >= 1")

    parts = partition_source_states(symmodel, nparts; strategy = partition_strategy)
    jobs = [
        (LocalGridBasedSymbolicModel(symmodel, Xset_local), concrete_system_approx) for
        Xset_local in parts
    ]

    results = if isempty(procs)
        map(jobs) do job
            local_symmodel, local_approx = job
            return _run_local_partition(
                local_symmodel,
                local_approx;
                verbose = false,
                update_interval = update_interval,
                progress_dt = progress_dt,
                threaded = threaded_per_worker,
            )
        end
    else
        pool = Distributed.WorkerPool(procs)
        Distributed.pmap(pool, jobs) do job
            local_symmodel, local_approx = job
            return _run_local_partition(
                local_symmodel,
                local_approx;
                verbose = false,
                update_interval = update_interval,
                progress_dt = progress_dt,
                threaded = threaded_per_worker,
            )
        end
    end

    transitions = Tuple{Int, Int, Int}[]
    total_sources = 0
    total_transitions = 0

    for res in results
        append!(transitions, res.transitions)
        total_sources += res.n_source_states
        total_transitions += res.n_transitions
    end

    verbose && @info(
        "Distributed abstraction finished",
        nparts = nparts,
        total_source_states = total_sources,
        total_transitions = total_transitions,
    )

    return transitions
end
