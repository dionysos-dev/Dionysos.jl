# ----------------------------
# AbstractStateSet interface
# ----------------------------

"""
    AbstractStateSet{N}

A set of abstract states (integer labels) over an `N`-dimensional concrete
space. The labels only have meaning relative to an [`AbstractMapping`](@ref),
which is why the membership/enumeration methods take the mapping `m` as an
argument. For a self-contained object that bundles the set with its mapping, use
[`MappedStateSet`](@ref) (the public surface); the `(S, m, …)` methods below are
the implementation it forwards to.

# Extending

Implement `contains_state`, `enum_states`, `add_state!`, `remove_state!`,
`empty_states!`. `get_n_state`, `add_states!`, `add_set!`, `remove_set!` have
generic defaults.
"""
abstract type AbstractStateSet{N} end

contains_state(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("implement `contains_state` for $(typeof(S))")
enum_states(S::AbstractStateSet{N}, m::AbstractMapping{N}) where {N} =
    error("implement `enum_states` for $(typeof(S))")
get_n_state(S::AbstractStateSet{N}, m::AbstractMapping{N}) where {N} =
    length(enum_states(S, m))

add_state!(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("implement `add_state!` for $(typeof(S))")
remove_state!(S::AbstractStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("implement `remove_state!` for $(typeof(S))")
empty_states!(S::AbstractStateSet{N}) where {N} =
    error("implement `empty_states!` for $(typeof(S))")

function add_states!(S::AbstractStateSet{N}, m::AbstractMapping{N}, states) where {N}
    for q in states
        add_state!(S, m, q)
    end
    return S
end

function add_set!(
    S::AbstractStateSet{N},
    m::AbstractMapping{N},
    set,
    incl_mode::INCL_MODE,
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    add_states!(S, m, states)
    return collect(states)
end

function stateset_from_states(
    ::Type{S},
    m::AbstractMapping{N},
    states,
) where {N, S <: AbstractStateSet{N}}
    out = S()
    add_states!(out, m, states)
    return out
end

# Convenient default: ExplicitIdSet (BitSet)
stateset_from_states(m::AbstractMapping{N}, states) where {N} =
    stateset_from_states(ExplicitIdSet{N}, m, states)

function remove_set!(
    S::AbstractStateSet{N},
    m::AbstractMapping{N},
    set,
    incl_mode::INCL_MODE,
) where {N}
    states = get_states_from_set(m, set, incl_mode)
    for q in states
        remove_state!(S, m, q)
    end
    return collect(states)
end
