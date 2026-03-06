import Base: intersect

# ----------------------------
# Core types / invariants
# ----------------------------

abstract type AbstractLazySet{N, T} end

get_dims(::AbstractLazySet{N, T}) where {N, T} = N

abstract type AbstractSetNode{N, T} <: AbstractLazySet{N, T} end

# ----------------------------
# LazySetUnion: contains ONLY nodes
# ----------------------------

struct LazySetUnion{N, T} <: AbstractLazySet{N, T}
    sets::Vector{AbstractSetNode{N, T}}
end

LazySetUnion{N, T}() where {N, T} = LazySetUnion{N, T}(AbstractSetNode{N, T}[])
LazySetUnion(v::Vector{<:AbstractSetNode{N, T}}) where {N, T} =
    LazySetUnion{N, T}(AbstractSetNode{N, T}[s for s in v])

Base.isempty(U::LazySetUnion) = isempty(U.sets)
get_sets(U::LazySetUnion) = U.sets

# ----------------------------
# LazySetIntersection:
# ----------------------------

struct LazySetIntersection{N, T} <: AbstractLazySet{N, T}
    sets::Vector{AbstractSetNode{N, T}}
end

LazySetIntersection{N, T}() where {N, T} =
    LazySetIntersection{N, T}(AbstractSetNode{N, T}[])

_flat_sets(s::AbstractLazySet{N, T}) where {N, T} = AbstractLazySet{N, T}[s]
_flat_sets(I::LazySetIntersection{N, T}) where {N, T} = I.sets

# ----------------------------
# LazySetMinus: children can be Node or LazySetUnion (but never LazySetMinus)
# ----------------------------

const MinusChild{N, T} = Union{AbstractSetNode{N, T}, LazySetUnion{N, T}}

struct LazySetMinus{N, T, CA <: MinusChild{N, T}, CB <: MinusChild{N, T}} <:
       AbstractLazySet{N, T}
    A::CA
    B::CB
end

# ----------------------------
# Guards / flatten helpers
# ----------------------------

_assert_not_minus(s) =
    s isa LazySetMinus ? throw(ArgumentError("Nested LazySetMinus forbidden")) : s

# Promote a child (node) to a union when we need to add more stuff to it
_promote_copy(child::AbstractSetNode{N, T}) where {N, T} = begin
    U = LazySetUnion{N, T}()
    push!(U.sets, child)
    U
end
_promote_copy(child::LazySetUnion{N, T}) where {N, T} = LazySetUnion{N, T}(copy(child.sets))

# Iterate nodes from either a node or a union
_eachnode(s::AbstractSetNode{N, T}) where {N, T} = (s,)
_eachnode(U::LazySetUnion{N, T}) where {N, T} = U.sets

# ----------------------------
# Membership
# ----------------------------

point_in_set(S, x) = (x ∈ S)

function point_in_set(U::LazySetUnion{N, T}, x) where {N, T}
    @inbounds for s in U.sets
        point_in_set(s, x) && return true
    end
    return false
end

point_in_set(S::LazySetMinus, x) = point_in_set(S.A, x) && !point_in_set(S.B, x)

Base.in(x, U::LazySetUnion) = point_in_set(U, x)
Base.in(x, S::LazySetMinus) = point_in_set(S, x)
Base.in(x, I::LazySetIntersection{N, T}) where {N, T} = all(s -> (x in s), I.sets)

# ----------------------------
# add_set! / remove_set!
# ----------------------------

function add_set!(U::LazySetUnion{N, T}, s) where {N, T}
    _assert_not_minus(s)
    if s isa AbstractSetNode{N, T}
        push!(U.sets, s)
    elseif s isa LazySetUnion{N, T}
        Base.append!(U.sets, s.sets)
    else
        throw(ArgumentError("Invalid type $(typeof(s)) for LazySetUnion{$N,$T}"))
    end
    return U
end

add_set(U::LazySetUnion{N, T}, s) where {N, T} =
    (V = LazySetUnion{N, T}(copy(U.sets)); add_set!(V, s); V)

function add_set(S::LazySetMinus{N, T}, s::LazySetMinus{N, T}) where {N, T}
    if isempty(s.B) || isequal(S.B, s.B)
        Aunion = _promote_copy(S.A)
        add_set!(Aunion, s.A)
        return LazySetMinus(Aunion, S.B)
    else
        throw(
            ArgumentError(
                "Nested LazySetMinus forbidden: cannot canonically add $(typeof(s)) to LazySetMinus unless s.B is empty or S.B == s.B",
            ),
        )
    end
end

function add_set(S::LazySetMinus{N, T}, s) where {N, T}
    _assert_not_minus(s)
    Aunion = _promote_copy(S.A)
    add_set!(Aunion, s)
    return LazySetMinus(Aunion, S.B)
end

function remove_set(S::LazySetMinus{N, T}, s) where {N, T}
    _assert_not_minus(s)
    Bunion = _promote_copy(S.B)
    add_set!(Bunion, s)
    return LazySetMinus(S.A, Bunion)
end

# ----------------------------
# Others
# ----------------------------

# LazySetUnion is a container of nodes
Base.length(U::LazySetUnion) = length(U.sets)
Base.eltype(::Type{LazySetUnion{N, T}}) where {N, T} = AbstractSetNode{N, T}
Base.iterate(U::LazySetUnion, st::Int = 1) =
    st > length(U.sets) ? nothing : (U.sets[st], st+1)

# Let LazySetMinus behave like "its A part" for iteration/length
Base.length(M::LazySetMinus{N, T}) where {N, T} = begin
    A = M.A
    A isa LazySetUnion{N, T} ? length(A) : 1
end
Base.eltype(::Type{LazySetMinus{N, T}}) where {N, T} = AbstractSetNode{N, T}
Base.iterate(M::LazySetMinus{N, T}, st::Int = 1) where {N, T} = begin
    A = M.A
    A isa LazySetUnion{N, T} ? iterate(A, st) : (st == 1 ? (A, 2) : nothing)
end

# ----------------------------
# Intersection
# ----------------------------

function _push_nonempty!(out::LazySetUnion{N, T}, s) where {N, T}
    _assert_not_minus(s)
    if s isa AbstractSetNode{N, T}
        isempty(s) || push!(out.sets, s)
    elseif s isa LazySetUnion{N, T}
        isempty(s) || append!(out.sets, s.sets)
    else
        throw(ArgumentError("Invalid type $(typeof(s)) for LazySetUnion{$N,$T}"))
    end
    return out
end

function intersect(U::LazySetUnion{N, T}, s::AbstractSetNode{N, T}) where {N, T}
    out = LazySetUnion{N, T}()
    @inbounds for a in U.sets
        _push_nonempty!(out, intersect(a, s))
    end
    return out
end

intersect(s::AbstractSetNode{N, T}, U::LazySetUnion{N, T}) where {N, T} = intersect(U, s)

# (A \ B) ∩ S = (A ∩ S) \ (B ∩ S)  [tighter than keeping B unchanged]
function intersect(M::LazySetMinus{N, T}, s::AbstractSetNode{N, T}) where {N, T}
    A2 = intersect(M.A, s)
    B2 = intersect(M.B, s)
    _assert_not_minus(A2);
    _assert_not_minus(B2)

    # If A2 is empty => empty
    isempty(A2) && return LazySetUnion{N, T}()

    return LazySetMinus(A2, B2)
end

intersect(s::AbstractSetNode{N, T}, M::LazySetMinus{N, T}) where {N, T} = intersect(M, s)

# U ∩ V = ⋃_{a∈U} ⋃_{b∈V} (a ∩ b)
function intersect(U::LazySetUnion{N, T}, V::LazySetUnion{N, T}) where {N, T}
    out = LazySetUnion{N, T}()
    # iterate over smaller outer loop to reduce overhead a bit
    A, B = length(U.sets) <= length(V.sets) ? (U, V) : (V, U)
    @inbounds for a in A.sets, b in B.sets
        _push_nonempty!(out, intersect(a, b))
    end
    return out
end

# U ∩ (A \ B) = (U ∩ A) \ (U ∩ B)
function intersect(U::LazySetUnion{N, T}, M::LazySetMinus{N, T}) where {N, T}
    A2 = intersect(U, M.A)  # returns LazySetUnion
    B2 = intersect(U, M.B)  # returns LazySetUnion
    isempty(A2) && return LazySetUnion{N, T}()
    return LazySetMinus(A2, B2)
end

# (A \ B) ∩ U = (A ∩ U) \ (B ∩ U)
function intersect(M::LazySetMinus{N, T}, U::LazySetUnion{N, T}) where {N, T}
    A2 = intersect(M.A, U)  # returns LazySetUnion
    B2 = intersect(M.B, U)  # returns LazySetUnion
    isempty(A2) && return LazySetUnion{N, T}()
    return LazySetMinus(A2, B2)
end

# (A \ B) ∩ (C \ D) = (A ∩ C) \ (B ∪ D)
# (safe over-approx of hole union; exact for set difference algebra)
function intersect(M1::LazySetMinus{N, T}, M2::LazySetMinus{N, T}) where {N, T}
    A2 = intersect(M1.A, M2.A)                 # allowed: child ∩ child
    isempty(A2) && return LazySetUnion{N, T}()
    Bunion = LazySetUnion{N, T}()                  # build B ∪ D
    add_set!(Bunion, M1.B)
    add_set!(Bunion, M2.B)
    return LazySetMinus(A2, Bunion)
end

# intersection ∩ anything
function intersect(I::LazySetIntersection{N, T}, s::AbstractLazySet{N, T}) where {N, T}
    return LazySetIntersection{N, T}(vcat(I.sets, _flat_sets(s)))
end
intersect(s::AbstractLazySet{N, T}, I::LazySetIntersection{N, T}) where {N, T} =
    intersect(I, s)

# ----------------------------
# Periodic wrapping
# ----------------------------

function set_in_period(U::LazySetUnion{N, T}, periodic_dims, periods, start) where {N, T}
    out = LazySetUnion{N, T}()
    for s in U.sets
        ws = set_in_period(s, periodic_dims, periods, start)
        _assert_not_minus(ws)
        add_set!(out, ws)
    end
    return out
end

function set_in_period(S::LazySetMinus{N, T}, periodic_dims, periods, start) where {N, T}
    A2 = set_in_period(S.A, periodic_dims, periods, start)
    B2 = set_in_period(S.B, periodic_dims, periods, start)
    _assert_not_minus(A2);
    _assert_not_minus(B2)
    return LazySetMinus(A2, B2)
end

# ----------------------------
# Plots
# ----------------------------

@recipe function f(U::LazySetUnion; dims = [1, 2], label = "set")
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

compare_sets(set1, set2) = get_volume(set1) < get_volume(set2)
function Base.sort!(
    iset::LazySetIntersection{N, T};
    lt = compare_sets,
    rev = false,
) where {N, T}
    sort!(iset.sets; lt = lt, rev = rev)
    return iset
end
@recipe function f(iset::LazySetIntersection)
    sets = sort(iset.sets; lt = compare_sets, rev = true)
    for set in sets
        @series begin
            return set
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
