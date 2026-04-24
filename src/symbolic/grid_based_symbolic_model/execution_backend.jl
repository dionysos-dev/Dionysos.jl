# ------------------------------------------------------------
# Execution Backends for Grid-Based Abstraction
# ------------------------------------------------------------

import Distributed

# ------------------------------------------------------------
# Execution Backend types
# ------------------------------------------------------------

"""
    AbstractExecutionBackend

Abstract type for execution backends used to compute the transition relation
in grid-based symbolic abstraction.

An execution backend defines **how the abstraction computation is executed**:
- sequentially,
- multithreaded,
- distributed across Julia workers,
- or via SLURM array jobs.
"""
abstract type AbstractExecutionBackend end

"""
    SequentialBackend()

Sequential execution (no parallelism).
"""
struct SequentialBackend <: AbstractExecutionBackend end

"""
    ThreadedBackend(progress_dt=0.2)

Multithreaded execution using all available Julia threads.

# Parameters
- `progress_dt`: minimum time (in seconds) between progress updates.
"""
struct ThreadedBackend <: AbstractExecutionBackend
    progress_dt::Float64
end

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

"""
    SlurmArrayBackend(
        nchunks,
        chunk_id=nothing,
        outdir,
        partition_strategy=:contiguous,
        write_only=true,
    )

Execution using SLURM array jobs (file-based parallelism).

# Parameters
- `nchunks`: total number of chunks (slurm array size).
- `chunk_id`: current chunk (defaults to `SLURM_ARRAY_TASK_ID`).
- `outdir`: directory where results are written.
- `partition_strategy`: how to split states (`:contiguous` or `:roundrobin`).
- `write_only`: if `true`, only writes transitions to disk (recommended for SLURM).
"""
struct SlurmArrayBackend <: AbstractExecutionBackend
    nchunks::Int
    chunk_id::Union{Nothing, Int}
    outdir::String
    partition_strategy::Symbol
    write_only::Bool
end

ThreadedBackend() = ThreadedBackend(0.2)

JuliaDistributedBackend() =
    JuliaDistributedBackend(nothing, nothing, :roundrobin, false, true)

SlurmArrayBackend(nchunks::Int, outdir::String) =
    SlurmArrayBackend(nchunks, nothing, outdir, :contiguous, true)

# ------------------------------------------------------------
# Local partition symbolic model
# ------------------------------------------------------------

"""
    LocalGridBasedSymbolicModel

Wrapper around a global symbolic model that overrides only the source domain.
The state mapping, input mapping, retained domain and input domain remain global.
"""
struct LocalGridBasedSymbolicModel{N, M, SM, XS} <: GridBasedSymbolicModel{N, M}
    parent::SM
    Xset_local::XS
end

function LocalGridBasedSymbolicModel(
    symmodel::GridBasedSymbolicModel{N, M},
    Xset_local,
) where {N, M}
    return LocalGridBasedSymbolicModel{N, M, typeof(symmodel), typeof(Xset_local)}(
        symmodel,
        Xset_local,
    )
end

get_state_mapping(sym::LocalGridBasedSymbolicModel) = get_state_mapping(sym.parent)

get_input_mapping(sym::LocalGridBasedSymbolicModel) = get_input_mapping(sym.parent)

get_source_domain(sym::LocalGridBasedSymbolicModel) = sym.Xset_local

get_retained_domain(sym::LocalGridBasedSymbolicModel) = get_retained_domain(sym.parent)

get_input_domain(sym::LocalGridBasedSymbolicModel) = get_input_domain(sym.parent)

get_concrete_state(sym::LocalGridBasedSymbolicModel, q) = get_concrete_state(sym.parent, q)

get_concrete_elem(sym::LocalGridBasedSymbolicModel, q::Int) =
    get_concrete_elem(sym.parent, q)

get_concrete_input(sym::LocalGridBasedSymbolicModel, u) = get_concrete_input(sym.parent, u)

get_abstract_state(sym::LocalGridBasedSymbolicModel, x) = get_abstract_state(sym.parent, x)

add_transitions!(sym::LocalGridBasedSymbolicModel, trans) =
    add_transitions!(sym.parent, trans)

# ------------------------------------------------------------
# Partitioning utilities
# ------------------------------------------------------------

function partition_source_state_ids(
    symmodel::GridBasedSymbolicModel,
    nparts::Int;
    strategy::Symbol = :roundrobin,
)
    nparts >= 1 || error("nparts must be >= 1")

    states = collect(enum_source_states(symmodel))
    parts = [Int[] for _ in 1:nparts]

    isempty(states) && return parts

    if strategy == :roundrobin
        for (k, q) in enumerate(states)
            push!(parts[mod1(k, nparts)], q)
        end

    elseif strategy == :contiguous
        n = length(states)
        for i in 1:nparts
            a = fld((i - 1) * n, nparts) + 1
            b = fld(i * n, nparts)
            a <= b && append!(parts[i], @view states[a:b])
        end

    else
        error("Unknown partition strategy: $strategy")
    end

    return parts
end

function local_symmodel_from_state_ids(
    symmodel::GridBasedSymbolicModel,
    state_ids::Vector{Int},
)
    Xmap = get_state_mapping(symmodel)
    Xset_local = MP.stateset_from_states(Xmap, state_ids)
    return LocalGridBasedSymbolicModel(symmodel, Xset_local)
end

# ------------------------------------------------------------
# Main execution backend dispatcher
# ------------------------------------------------------------

function compute_abstract_system_from_concrete_system!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    execution_backend::AbstractExecutionBackend = SequentialBackend(),
    kwargs...,
)
    return compute_abstract_system_from_concrete_system!(
        symmodel,
        concrete_system_approx,
        execution_backend;
        kwargs...,
    )
end

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

function compute_abstract_system_from_concrete_system!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx,
    execution_backend::SlurmArrayBackend;
    kwargs...,
)
    return compute_abstract_system_slurm_array!(
        symmodel,
        concrete_system_approx,
        execution_backend;
        kwargs...,
    )
end
