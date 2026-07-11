"""
    MappedStateSet{N, S, M}

A state set bundled with the mapping that gives its integer labels meaning.
Consumers hold one self-contained object instead of threading the
`(set, mapping)` pair through every call — which also removes the
"wrong mapping" bug class. This is the public surface for state sets; the
`(set, mapping, …)` methods on [`AbstractStateSet`](@ref) are what it forwards to.
"""
struct MappedStateSet{N, S <: AbstractStateSet{N}, M <: AbstractMapping{N}}
    set::S
    mapping::M
end

get_set(ms::MappedStateSet) = ms.set
get_mapping(ms::MappedStateSet) = ms.mapping
get_dim(::MappedStateSet{N}) where {N} = N

contains_state(ms::MappedStateSet, q::Int) = contains_state(ms.set, ms.mapping, q)
enum_states(ms::MappedStateSet) = enum_states(ms.set, ms.mapping)
get_n_state(ms::MappedStateSet) = get_n_state(ms.set, ms.mapping)

add_state!(ms::MappedStateSet, q::Int) = add_state!(ms.set, ms.mapping, q)
remove_state!(ms::MappedStateSet, q::Int) = remove_state!(ms.set, ms.mapping, q)
add_states!(ms::MappedStateSet, states) = add_states!(ms.set, ms.mapping, states)
add_set!(ms::MappedStateSet, set, incl_mode::INCL_MODE) =
    add_set!(ms.set, ms.mapping, set, incl_mode)
remove_set!(ms::MappedStateSet, set, incl_mode::INCL_MODE) =
    remove_set!(ms.set, ms.mapping, set, incl_mode)

Base.copy(ms::MappedStateSet) = MappedStateSet(copy(ms.set), ms.mapping)

@recipe function f(ms::MappedStateSet)
    return ((ms.set, ms.mapping),) # delegates to the tuple recipe
end
