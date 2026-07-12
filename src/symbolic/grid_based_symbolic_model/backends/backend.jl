# ------------------------------------------------------------
# Execution Backends for Grid-Based Abstraction
# ------------------------------------------------------------

import Distributed

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
