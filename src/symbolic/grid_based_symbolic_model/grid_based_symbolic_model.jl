"""
    GridBasedSymbolicModel{N,M} <: SymbolicModel{N,M}

Intermediate abstract type for symbolic models that rely on a grid-based or
mapping-based discretization.

Semantics:
- state mapping: global abstract-state numbering / coordinate map
- input mapping: global abstract-input numbering / coordinate map
- source domain (`Xset`): states enumerated as sources
- retained domain (`Rset`): states allowed as targets
"""
abstract type GridBasedSymbolicModel{N,M} <: SymbolicModel{N,M} end

# ------------------------------------------------------------------
# Required interface for concrete subtypes
# ------------------------------------------------------------------

get_source_domain(symmodel::GridBasedSymbolicModel) =
    error("get_source_domain not implemented for $(typeof(symmodel))")

get_retained_domain(symmodel::GridBasedSymbolicModel) =
    error("get_retained_domain not implemented for $(typeof(symmodel))")

# ------------------------------------------------------------------
# Default generic behavior
# ------------------------------------------------------------------

enum_source_states(symmodel::GridBasedSymbolicModel) =
    MP.enum_states(get_source_domain(symmodel), get_state_mapping(symmodel))

is_allowed_state(symmodel::GridBasedSymbolicModel, q) =
    q !== nothing &&
    MP.contains_state(get_retained_domain(symmodel), get_state_mapping(symmodel), q)

get_n_state(symmodel::GridBasedSymbolicModel) =
    length(collect(enum_source_states(symmodel)))

get_n_input(symmodel::GridBasedSymbolicModel) =
    length(collect(enum_inputs(symmodel)))

get_state_domain(symmodel::GridBasedSymbolicModel) = get_source_domain(symmodel)
enum_states(symmodel::GridBasedSymbolicModel) = enum_source_states(symmodel)


include("sequential_threaded_abstraction.jl")
include("distributed_abstraction.jl")

function compute_abstract_system_from_concrete_system!(
    symmodel::GridBasedSymbolicModel,
    concrete_system_approx;
    distributed::Bool = false,
    threaded::Bool = false,
    procs = Distributed.workers(),
    nparts::Int = max(length(procs), 1),
    partition_strategy::Symbol = :roundrobin,
    kwargs...,
)
    if distributed
        return compute_abstract_system_distributed!(
            symmodel,
            concrete_system_approx;
            procs = procs,
            nparts = nparts,
            partition_strategy = partition_strategy,
            threaded_per_worker = threaded,
            kwargs...,
        )
    else
        compute_abstract_system!(
            symmodel,
            concrete_system_approx;
            threaded = threaded,
            kwargs...,
        )
        return symmodel
    end
end
