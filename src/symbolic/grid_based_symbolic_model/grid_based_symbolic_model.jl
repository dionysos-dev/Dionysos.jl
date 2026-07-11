"""
    GridBasedSymbolicModel{N,M} <: SymbolicModel{N,M}

Intermediate abstract type for symbolic models that rely on a grid-based or
mapping-based discretization.

Semantics:
- state mapping: global abstract-state numbering / coordinate map
- input mapping: global abstract-input numbering / coordinate map
- state set (`Xset`): states enumerated as sources
- retained set (`Rset`): states allowed as targets
"""
abstract type GridBasedSymbolicModel{N, M} <: SymbolicModel{N, M} end

# ------------------------------------------------------------------
# Required interface for concrete subtypes
# ------------------------------------------------------------------

get_state_set(symmodel::GridBasedSymbolicModel) =
    error("get_state_set not implemented for $(typeof(symmodel))")

get_retained_set(symmodel::GridBasedSymbolicModel) =
    error("get_retained_set not implemented for $(typeof(symmodel))")

# ------------------------------------------------------------------
# Default generic behavior
# ------------------------------------------------------------------

# Intent-naming alias: the states enumerated as transition *sources* (as opposed
# to the retained set of allowed *targets*). Same as `enum_states`; used at the
# partitioning / kernel call sites where "source" carries meaning.
enum_source_states(symmodel::GridBasedSymbolicModel) = enum_states(symmodel)

is_allowed_state(symmodel::GridBasedSymbolicModel, q::Int) =
    q !== nothing && MP.contains_state(get_mapped_retained_set(symmodel), q)

@inline function _keep_state_input(
    symmodel::GridBasedSymbolicModel,
    abstract_state::Int,
    abstract_input::Int,
    state_input_filter::Union{Nothing, Function},
)
    state_input_filter === nothing && return true
    x = get_concrete_state(symmodel, abstract_state)
    u = get_concrete_input(symmodel, abstract_input)
    return state_input_filter(x, u)
end

function filtered_source_states(
    symmodel::GridBasedSymbolicModel,
    state_filter::Union{Nothing, Function},
)
    states = collect(enum_source_states(symmodel))
    state_filter === nothing && return states
    return [q for q in states if state_filter(get_concrete_state(symmodel, q))]
end
