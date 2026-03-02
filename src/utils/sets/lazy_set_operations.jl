import Base: intersect

# ----------------------------
# Core types / invariants
# ----------------------------

abstract type AbstractLazySet{N,T} end

get_dims(::AbstractLazySet{N,T}) where {N,T} = N

abstract type AbstractSetNode{N,T} <: AbstractLazySet{N,T} end

# ----------------------------
# LazyUnion: contains ONLY nodes
# ----------------------------

struct LazyUnion{N,T} <: AbstractLazySet{N,T}
    sets::Vector{AbstractSetNode{N,T}}
end

LazyUnion{N,T}() where {N,T} = LazyUnion{N,T}(AbstractSetNode{N,T}[])
LazyUnion(v::Vector{<:AbstractSetNode{N,T}}) where {N,T} =
    LazyUnion{N,T}(AbstractSetNode{N,T}[s for s in v])

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
        Base.append!(U.sets, s.sets)
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
# Others
# ----------------------------

# LazyUnion is a container of nodes
Base.length(U::LazyUnion) = length(U.sets)
Base.eltype(::Type{LazyUnion{N,T}}) where {N,T} = AbstractSetNode{N,T}
Base.iterate(U::LazyUnion, st::Int=1) = st > length(U.sets) ? nothing : (U.sets[st], st+1)

# Let LazySetMinus behave like "its A part" for iteration/length
Base.length(M::LazySetMinus{N,T}) where {N,T} = begin
    A = M.A
    A isa LazyUnion{N,T} ? length(A) : 1
end
Base.eltype(::Type{LazySetMinus{N,T}}) where {N,T} = AbstractSetNode{N,T}
Base.iterate(M::LazySetMinus{N,T}, st::Int=1) where {N,T} = begin
    A = M.A
    A isa LazyUnion{N,T} ? iterate(A, st) : (st == 1 ? (A, 2) : nothing)
end

# ----------------------------
# Intersection
# ----------------------------

function _push_nonempty!(out::LazyUnion{N,T}, s) where {N,T}
    _assert_not_minus(s)
    if s isa AbstractSetNode{N,T}
        isempty(s) || push!(out.sets, s)
    elseif s isa LazyUnion{N,T}
        isempty(s) || append!(out.sets, s.sets)
    else
        throw(ArgumentError("Invalid type $(typeof(s)) for LazyUnion{$N,$T}"))
    end
    return out
end

function intersect(U::LazyUnion{N,T}, s::AbstractSetNode{N,T}) where {N,T}
    out = LazyUnion{N,T}()
    @inbounds for a in U.sets
        _push_nonempty!(out, intersect(a, s))
    end
    return out
end

intersect(s::AbstractSetNode{N,T}, U::LazyUnion{N,T}) where {N,T} = intersect(U, s)

# (A \ B) ∩ S = (A ∩ S) \ (B ∩ S)  [tighter than keeping B unchanged]
function intersect(M::LazySetMinus{N,T}, s::AbstractSetNode{N,T}) where {N,T}
    A2 = intersect(M.A, s)
    B2 = intersect(M.B, s)
    _assert_not_minus(A2); _assert_not_minus(B2)

    # If A2 is empty => empty
    isempty(A2) && return LazyUnion{N,T}()

    return LazySetMinus(A2, B2)
end

intersect(s::AbstractSetNode{N,T}, M::LazySetMinus{N,T}) where {N,T} = intersect(M, s)

# U ∩ V = ⋃_{a∈U} ⋃_{b∈V} (a ∩ b)
function intersect(U::LazyUnion{N,T}, V::LazyUnion{N,T}) where {N,T}
    out = LazyUnion{N,T}()
    # iterate over smaller outer loop to reduce overhead a bit
    A, B = length(U.sets) <= length(V.sets) ? (U, V) : (V, U)
    @inbounds for a in A.sets, b in B.sets
        _push_nonempty!(out, intersect(a, b))
    end
    return out
end

# U ∩ (A \ B) = (U ∩ A) \ (U ∩ B)
function intersect(U::LazyUnion{N,T}, M::LazySetMinus{N,T}) where {N,T}
    A2 = intersect(U, M.A)  # returns LazyUnion
    B2 = intersect(U, M.B)  # returns LazyUnion
    isempty(A2) && return LazyUnion{N,T}()
    return LazySetMinus(A2, B2)
end

# (A \ B) ∩ U = (A ∩ U) \ (B ∩ U)
function intersect(M::LazySetMinus{N,T}, U::LazyUnion{N,T}) where {N,T}
    A2 = intersect(M.A, U)  # returns LazyUnion
    B2 = intersect(M.B, U)  # returns LazyUnion
    isempty(A2) && return LazyUnion{N,T}()
    return LazySetMinus(A2, B2)
end

# (A \ B) ∩ (C \ D) = (A ∩ C) \ (B ∪ D)
# (safe over-approx of hole union; exact for set difference algebra)
function intersect(M1::LazySetMinus{N,T}, M2::LazySetMinus{N,T}) where {N,T}
    A2 = intersect(M1.A, M2.A)                 # allowed: child ∩ child
    isempty(A2) && return LazyUnion{N,T}()
    Bunion = LazyUnion{N,T}()                  # build B ∪ D
    add_set!(Bunion, M1.B)
    add_set!(Bunion, M2.B)
    return LazySetMinus(A2, Bunion)
end

# ----------------------------
# Periodic wrapping
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

# ----------------------------
# Plots
# ----------------------------

@recipe function f(U::LazyUnion; dims = [1,2], label = "set")
    dims := dims

    first_series = true
    for s in U.sets
        @series begin
            label := first_series ? label : ""
            first_series = false
            return s
        end
    end
end

@recipe function f(
    S::LazySetMinus;
    dims = [1, 2],
    hole_color = :gray,
    hole_alpha = 1.0,
    label = "Set",
)
    dims := dims

    @series begin
        label := label
        return S.A
    end

    @series begin
        label := ""
        seriestype := :shape
        fillcolor := hole_color
        fillalpha := hole_alpha
        return S.B
    end
end