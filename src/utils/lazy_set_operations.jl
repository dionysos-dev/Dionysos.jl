# ----------------------------
# Core types / invariants
# ----------------------------

abstract type AbstractLazySet{N,T} end
abstract type AbstractSetNode{N,T} <: AbstractLazySet{N,T} end

# ----------------------------
# LazyUnion: contains ONLY nodes
# ----------------------------

struct LazyUnion{N,T} <: AbstractLazySet{N,T}
    sets::Vector{AbstractSetNode{N,T}}
end

LazyUnion{N,T}() where {N,T} = LazyUnion{N,T}(AbstractSetNode{N,T}[])

Base.isempty(U::LazyUnion) = isempty(U.sets)
get_sets(U::LazyUnion) = U.sets

# ----------------------------
# LazySetMinus: children can be Node or LazyUnion (but never LazySetMinus)
# ----------------------------

const MinusChild{N,T} = Union{AbstractSetNode{N,T}, LazyUnion{N,T}}

struct LazySetMinus{N,T,CA<:MinusChild{N,T},CB<:MinusChild{N,T}} <: AbstractLazySet{N,T}
    A::CA
    B::CB
end

LazySetMinus(A::MinusChild{N,T}, B::MinusChild{N,T}) where {N,T} =
    LazySetMinus{N,T,typeof(A),typeof(B)}(A, B)

# ----------------------------
# Guards / flatten helpers
# ----------------------------

_assert_not_minus(s) = s isa LazySetMinus ? throw(ArgumentError("Nested LazySetMinus forbidden")) : s

# Promote a child (node) to a union when we need to add more stuff to it
_promote_copy(child::AbstractSetNode{N,T}) where {N,T} = begin
    U = LazyUnion{N,T}()
    push!(U.sets, child)
    U
end
_promote_copy(child::LazyUnion{N,T}) where {N,T} =
    LazyUnion{N,T}(copy(child.sets))

# Iterate nodes from either a node or a union
_eachnode(s::AbstractSetNode{N,T}) where {N,T} = (s,)
_eachnode(U::LazyUnion{N,T}) where {N,T} = U.sets

# ----------------------------
# Membership
# ----------------------------

point_in_set(S, x) = (x ∈ S)

function point_in_set(U::LazyUnion{N,T}, x) where {N,T}
    @inbounds for s in U.sets
        point_in_set(s, x) && return true
    end
    return false
end

point_in_set(S::LazySetMinus, x) =
    point_in_set(S.A, x) && !point_in_set(S.B, x)

Base.in(x, U::LazyUnion) = point_in_set(U, x)
Base.in(x, S::LazySetMinus) = point_in_set(S, x)

# ----------------------------
# add_set! / remove_set!
# ----------------------------

"""
    add_set!(U::LazyUnion, s)

Adds node(s) into union. Flattens unions. Rejects minus.
"""
function add_set!(U::LazyUnion{N,T}, s) where {N,T}
    _assert_not_minus(s)
    if s isa AbstractSetNode{N,T}
        push!(U.sets, s)
    elseif s isa LazyUnion{N,T}
        append!(U.sets, s.sets)
    else
        throw(ArgumentError("Invalid type $(typeof(s)) for LazyUnion{$N,$T}"))
    end
    return U
end

add_set(U::LazyUnion{N,T}, s) where {N,T} = (V = LazyUnion{N,T}(copy(U.sets)); add_set!(V,s); V)

function add_set(S::LazySetMinus{N,T}, s) where {N,T}
    _assert_not_minus(s)
    Aunion = _promote_copy(S.A)
    add_set!(Aunion, s)
    return LazySetMinus(Aunion, S.B)
end

function remove_set(S::LazySetMinus{N,T}, s) where {N,T}
    _assert_not_minus(s)
    Bunion = _promote_copy(S.B)
    add_set!(Bunion, s)
    return LazySetMinus(S.A, Bunion)
end

# ----------------------------
# Periodic wrapping (delegates to leaves)
# ----------------------------

function set_in_period(U::LazyUnion{N,T}, periodic_dims, periods, start) where {N,T}
    out = LazyUnion{N,T}()
    for s in U.sets
        ws = set_in_period(s, periodic_dims, periods, start)
        _assert_not_minus(ws)
        add_set!(out, ws)
    end
    return out
end

function set_in_period(S::LazySetMinus{N,T}, periodic_dims, periods, start) where {N,T}
    A2 = set_in_period(S.A, periodic_dims, periods, start)
    B2 = set_in_period(S.B, periodic_dims, periods, start)
    _assert_not_minus(A2); _assert_not_minus(B2)
    return LazySetMinus(A2, B2)
end