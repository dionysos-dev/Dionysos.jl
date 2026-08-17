# ============================================================
# Set algebra, backed by LazySets.
#
# All Dionysos sets ARE LazySets: boxes are `Hyperrectangle`, ellipsoids are
# `Ellipsoid`, and composite regions reuse the ecosystem's lazy wrappers:
#
#     union   ⋃ᵢ Aᵢ   → `LazySets.UnionSetArray`
#     minus   A \ B    → `LazySets.Intersection(A, LazySets.Complement(B))`
#
# On top we add only what LazySets lacks and Dionysos needs: exact ellipsoid
# kernels (`is_included`/`is_disjoint` specializations), periodic-domain
# wrapping (`set_in_period`) and incremental union/hole builders
# (`add_set`/`remove_set`, used to grow `ImplicitStateSet`s).
# ============================================================

"""
    is_included(X, Y) -> Bool

Whether `X ⊆ Y`. Generic sets delegate to LazySets; an exact analytic kernel
runs for two `LazySets.Ellipsoid`s (see `sets/ellipsoid_inclusion.jl`) —
Dionysos owns this verb because extending `Base.issubset` on two
LazySets-owned types would be piracy.
"""
is_included(X, Y) = X ⊆ Y

"""
    is_disjoint(X, Y) -> Bool

Whether `X ∩ Y = ∅`. Generic sets delegate to LazySets; an exact analytic
kernel runs for two `LazySets.Ellipsoid`s (see `sets/ellipsoid_intersection.jl`).
"""
is_disjoint(X, Y) = LazySets.isdisjoint(X, Y)

# ------------------------------------------------------------
# Axis-aligned boxes: `LazySets.Hyperrectangle` is the box type; an empty
# region is `LazySets.EmptySet` (never a crossed-bounds box).
# ------------------------------------------------------------

# Internal. `Hyperrectangle` pinned to `SVector` storage — note the parameter order is not
# `Hyperrectangle`'s. It exists so a container of boxes has a *concrete* element type, which
# propagates into `UnionSetArray`'s parameters and keeps membership statically dispatched
# (measured ~5× on `x ∈ U`). Never annotate an argument or field with it: it is one backing
# among several, and `LazySets.AbstractHyperrectangle` is what code should accept.
const _Box{N, T} = LazySets.Hyperrectangle{T, SVector{N, T}, SVector{N, T}}

# Full side lengths (2·radius); LazySets only exposes the half-widths.
get_h(H::LazySets.AbstractHyperrectangle) = 2 .* LazySets.radius_hyperrectangle(H)

"Sample a point of `X` uniformly (LazySets rejection sampling)."
sample(X::LazySets.LazySet) = SVector{LazySets.dim(X)}(LazySets.sample(X))

"Sample `N` points of `X` uniformly (LazySets rejection sampling); pass `rng` for reproducibility."
function samples(X::LazySets.LazySet, N::Int; rng = nothing)
    n = LazySets.dim(X)
    pts = rng === nothing ? LazySets.sample(X, N) : LazySets.sample(X, N; rng = rng)
    return [SVector{n}(p) for p in pts]
end

# Ellipsoids are sampled DIRECTLY — uniform in the unit ball, mapped through a
# square root of the shape matrix. Rejection sampling degenerates on the
# needle-thin ellipsoids funnel chains produce (LazySets' sampler then dies in
# a `vertices_list` fallback with no Ellipsoid method); the direct map is exact
# for any positive semi-definite shape, degenerate included.
function _sample_in_ellipsoid(E::LazySets.Ellipsoid, rng)
    c = collect(Float64, LazySets.center(E))
    n = length(c)
    F = eigen(Symmetric(Matrix{Float64}(LazySets.shape_matrix(E))))
    A = F.vectors * Diagonal(sqrt.(max.(F.values, 0.0)))
    u = rng === nothing ? randn(n) : randn(rng, n)
    nu = norm(u)
    nu == 0.0 && return SVector{n}(c)
    r = (rng === nothing ? rand() : rand(rng))^(1 / n)
    return SVector{n}(c .+ A * ((r / nu) .* u))
end

sample(E::LazySets.Ellipsoid) = _sample_in_ellipsoid(E, nothing)
samples(E::LazySets.Ellipsoid, N::Int; rng = nothing) =
    [_sample_in_ellipsoid(E, rng) for _ in 1:N]

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

_empty_like(r::LazySets.LazySet{T}) where {T} = LazySets.EmptySet{T}(LazySets.dim(r))
_empty_like(r) = LazySets.EmptySet(LazySets.dim(r))

# A `set_minus` is an `Intersection` whose second operand is a `Complement`;
# the pattern has no LazySets name, so Dionysos aliases it for dispatch.
const SetMinus{T, XT, BT} = LazySets.Intersection{T, XT, LazySets.Complement{T, BT}}

_is_empty_region(::LazySets.EmptySet) = true
_is_empty_region(B::LazySets.UnionSetArray) = isempty(B.array)
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
minus_hole(S) = _empty_like(S)

# ------------------------------------------------------------
# Incremental builders — grow the included region (`add_set`) or the holes
# (`remove_set`) of a `set_minus`. Used by `ImplicitStateSet`.
# ------------------------------------------------------------

_num_type(::LazySets.LazySet{T}) where {T} = T

# `UnionSetArray` needs a `Vector{<:LazySet{T}}`, not the unparametrised `LazySet`.
_region_vector(::Type{T}) where {T} = LazySets.LazySet{T}[]

_collect_sets!(arr, u::LazySets.UnionSetArray) = (append!(arr, u.array); arr)
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

# Tight outer bounding box; exact for boxes (identity) and any set with a
# support function (via LazySets), with cheap structural cases for composites.
_outer_box(S::LazySets.LazySet) = LazySets.box_approximation(S)
_outer_box(H::LazySets.AbstractHyperrectangle) = H

function _outer_box(U::LazySets.UnionSetArray)
    isempty(U.array) && error("cannot compute outer box of an empty union")
    boxes = [_outer_box(s) for s in U.array]
    lb = reduce((a, b) -> min.(a, b), (LazySets.low(B) for B in boxes))
    ub = reduce((a, b) -> max.(a, b), (LazySets.high(B) for B in boxes))
    return LazySets.Hyperrectangle(; low = lb, high = ub)
end
_outer_box(S::LazySets.Intersection) = _outer_box(minus_included(S))

# ------------------------------------------------------------
# Periodic wrapping (splits sets crossing a period boundary; no LazySets
# equivalent).
# ------------------------------------------------------------

# Wrap the interval [lb, ub] into [T0, T0 + Tper); returns one interval, or two
# when the wrapped interval straddles the period boundary.
function one_direction(lb, ub, Tper, T0)
    if ub - lb >= Tper
        return [(T0, T0 + Tper)]
    else
        lbw = T0 + mod(lb - T0, Tper)
        ubw = T0 + mod(ub - T0, Tper)
        if lbw <= ubw
            return [(lbw, ubw)]
        else
            return [(T0, ubw), (lbw, T0 + Tper)]
        end
    end
end

function _recursive_period_split!(
    out,
    rlow::NTuple{N},
    rhigh::NTuple{N},
    lb::NTuple{N},
    ub::NTuple{N},
    periodic_dims::SVector{P, Int},
    periods::SVector{P},
    start::SVector{P},
    i::Int,
) where {N, P}
    if i > P
        push!(out, LazySets.Hyperrectangle(; low = SVector(lb), high = SVector(ub)))
        return
    end
    dim = periodic_dims[i]
    intervals = one_direction(rlow[dim], rhigh[dim], periods[i], start[i])
    for (lI, uI) in intervals
        l = ntuple(k -> k == dim ? lI : lb[k], Val(N))
        u = ntuple(k -> k == dim ? uI : ub[k], Val(N))
        _recursive_period_split!(
            out,
            rlow,
            rhigh,
            l,
            u,
            periodic_dims,
            periods,
            start,
            i + 1,
        )
    end
    return
end

"""
    set_in_period(rect, periodic_dims, periods, start) -> LazySets.UnionSetArray

Split `rect` along periodic boundaries and return the union of the wrapped boxes.
"""
function set_in_period(
    rect::LazySets.AbstractHyperrectangle,
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T},
    start::SVector{P, T},
) where {P, T}
    N = LazySets.dim(rect)
    lb = ntuple(i -> LazySets.low(rect, i), N)
    ub = ntuple(i -> LazySets.high(rect, i), N)
    L = _Box{N, T}[]
    _recursive_period_split!(L, lb, ub, lb, ub, periodic_dims, periods, start, 1)
    return set_union(L)
end

function set_in_period(U::LazySets.UnionSetArray, periodic_dims, periods, start)
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

"Project `S` onto the coordinates `dims` (identity when `dims` covers `S`)."
project_set(S, dims) =
    collect(dims) == collect(1:LazySets.dim(S)) ? S : LazySets.project(S, collect(dims))

@recipe function f(U::LazySets.UnionSetArray; dims = [1, 2], label = "set")
    first_series = true
    for s in U.array
        @series begin
            label := first_series ? label : ""
            first_series = false
            return project_set(s, dims)
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
    @series begin
        label := label
        return project_set(minus_included(S), dims)
    end
    @series begin
        label := ""
        seriestype := :shape
        fillcolor := hole_color
        fillalpha := hole_alpha
        return project_set(minus_hole(S), dims)
    end
end
