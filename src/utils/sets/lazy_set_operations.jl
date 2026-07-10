# ============================================================
# Set algebra, backed by LazySets.
#
# Dionysos leaf sets (rectangles, ellipsoids, …) are `LazySets.LazySet`s (see
# `AbstractSetNode`). Composite regions reuse the ecosystem's lazy wrappers:
#
#     union   ⋃ᵢ Aᵢ   → `LazySets.UnionSetArray`
#     minus   A \ B    → `LazySets.Intersection(A, LazySets.Complement(B))`
#
# On top we add only what LazySets lacks and Dionysos needs: periodic-domain
# wrapping (`set_in_period`) and incremental union/hole builders
# (`add_set`/`remove_set`, used to grow `ImplicitStateSet`s).
# ============================================================

# The Dionysos leaf-set hierarchy is rooted in `LazySets.LazySet` so leaves
# compose with the ecosystem's set algebra. `N` = ambient dimension, `T` = the
# coordinate type (the `LazySet` numeric parameter).
abstract type AbstractLazySet{N, T} <: LazySets.LazySet{T} end

get_dim(::AbstractLazySet{N, T}) where {N, T} = N
get_dim(s::LazySets.LazySet) = LazySets.dim(s)
LazySets.dim(::AbstractLazySet{N, T}) where {N, T} = N

abstract type AbstractSetNode{N, T} <: AbstractLazySet{N, T} end
_outer_box(X::AbstractSetNode) = error("implement `_outer_box` for $(typeof(X))")

# ------------------------------------------------------------
# Composite-set constructors and dispatch aliases
# ------------------------------------------------------------

"""
    set_union(sets) -> LazySets.UnionSetArray

Lazy union `⋃ᵢ setsᵢ` of the given (Dionysos or LazySets) sets.
"""
set_union(sets::AbstractVector) = LazySets.UnionSetArray(collect(sets))
set_union(sets...) = set_union(collect(sets))

"""
    set_minus(A, B) -> LazySets.Intersection

The set difference `A \\ B`, represented lazily as `A ∩ Bᶜ`.
"""
function set_minus(A, B)
    # `A \ ∅ = A`; also avoids building an `Intersection` with a dimensionless
    # empty union (an empty `UnionSetArray` has no inferable dimension).
    _is_empty_region(B) && return A
    return LazySets.Intersection(A, LazySets.Complement(B))
end

"""
    empty_region(::Type{T}, n) -> LazySets.EmptySet{T}

An empty set of dimension `n` (a typed identity for `set_union`/`set_minus`; an
empty `UnionSetArray` has no well-defined dimension, so we use `EmptySet`).
"""
empty_region(::Type{T}, n::Integer) where {T} = LazySets.EmptySet{T}(n)
empty_region(n::Integer) = empty_region(Float64, n)
empty_region_like(r::AbstractLazySet{N, T}) where {N, T} = empty_region(T, N)
empty_region_like(r) = empty_region(Float64, LazySets.dim(r))

set_dim(s) = LazySets.dim(s)

# A `set_union` is any `UnionSetArray`; a `set_minus` is an `Intersection` whose
# second operand is a `Complement`.
const SetUnion = LazySets.UnionSetArray
const SetMinus{T, XT, BT} = LazySets.Intersection{T, XT, LazySets.Complement{T, BT}}
const EmptyRegion = LazySets.EmptySet
const Region = LazySets.LazySet

_is_empty_region(::LazySets.EmptySet) = true
_is_empty_region(B::SetUnion) = isempty(B.array)
_is_empty_region(B) = false

# LazySets' smart constructors may simplify `set_minus` away from an `Intersection`
# (e.g. `∅ \ B → ∅`, `A \ ∅ → A`), so these extractors are total: a plain set is
# "the included region with no hole".
"Included region `A` of `A \\ B`."
minus_included(S::LazySets.Intersection) = S.X
minus_included(S::LazySets.EmptySet) = S
minus_included(S) = S
"Hole `B` of `A \\ B`."
minus_hole(S::LazySets.Intersection) = S.Y.X
minus_hole(S::LazySets.EmptySet) = S
minus_hole(S) = empty_region_like(S)

# ------------------------------------------------------------
# Membership (LazySets provides `∈` for all of these; thin alias kept for the
# discretisation call sites that spell it out).
# ------------------------------------------------------------

point_in_set(S, x) = x ∈ S

# ------------------------------------------------------------
# Incremental builders — grow the included region (`add_set`) or the holes
# (`remove_set`) of a `set_minus`. Used by `ImplicitStateSet`.
# ------------------------------------------------------------

_num_type(::LazySets.LazySet{T}) where {T} = T

# `UnionSetArray` needs a `Vector{<:LazySet{T}}`, not the unparametrised `LazySet`.
_region_vector(::Type{T}) where {T} = LazySets.LazySet{T}[]

_collect_sets!(arr, u::SetUnion) = (append!(arr, u.array); arr)
_collect_sets!(arr, ::LazySets.EmptySet) = arr
_collect_sets!(arr, s) = (push!(arr, s); arr)

function _union_with(existing, s)
    arr = _region_vector(_num_type(s))
    _collect_sets!(arr, existing)
    _collect_sets!(arr, s)
    return LazySets.UnionSetArray(arr)
end

add_set(S, s) = set_minus(_union_with(minus_included(S), s), minus_hole(S))
remove_set(S, s) = set_minus(minus_included(S), _union_with(minus_hole(S), s))

# ------------------------------------------------------------
# Outer bounding box
# ------------------------------------------------------------

function _outer_box(U::SetUnion)
    isempty(U.array) && error("cannot compute outer box of an empty union")
    boxes = [_outer_box(s) for s in U.array]
    lb = reduce((a, b) -> min.(a, b), (B.lb for B in boxes))
    ub = reduce((a, b) -> max.(a, b), (B.ub for B in boxes))
    return HyperRectangle(lb, ub)
end
_outer_box(S::LazySets.Intersection) = _outer_box(minus_included(S))

# ------------------------------------------------------------
# Periodic wrapping (splits sets crossing a period boundary; no LazySets
# equivalent). `set_in_period` on a leaf `HyperRectangle` lives in rectangle.jl.
# ------------------------------------------------------------

function set_in_period(U::SetUnion, periodic_dims, periods, start)
    arr = _region_vector(_num_type(U))
    for s in U.array
        _collect_sets!(arr, set_in_period(s, periodic_dims, periods, start))
    end
    return LazySets.UnionSetArray(arr)
end

function set_in_period(S::LazySets.Intersection, periodic_dims, periods, start)
    A = set_in_period(minus_included(S), periodic_dims, periods, start)
    B = set_in_period(minus_hole(S), periodic_dims, periods, start)
    return set_minus(A, B)
end

set_in_period(S::LazySets.EmptySet, periodic_dims, periods, start) = S

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

@recipe function f(U::SetUnion; dims = [1, 2], label = "set")
    dims := dims
    first_series = true
    for s in U.array
        @series begin
            label := first_series ? label : ""
            first_series = false
            return s
        end
    end
end

@recipe function f(
    S::LazySets.Intersection{N, S1, S2};
    dims = [1, 2],
    hole_color = :gray,
    hole_alpha = 1.0,
    label = "Set",
) where {N, S1 <: LazySets.LazySet{N}, S2 <: LazySets.Complement{N}}
    dims := dims
    @series begin
        label := label
        return minus_included(S)
    end
    @series begin
        label := ""
        seriestype := :shape
        fillcolor := hole_color
        fillalpha := hole_alpha
        return minus_hole(S)
    end
end
