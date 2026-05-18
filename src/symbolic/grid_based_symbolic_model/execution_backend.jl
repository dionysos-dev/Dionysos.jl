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
