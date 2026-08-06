"""
    abstract type Grid{N, T} end

Defines an abstract type for grid-based structures in `N` dimensions with floating-point values `T`.
"""
abstract type Grid{N, T} end

# ----------------------------
# Required API for Grids
# ----------------------------
get_origin(grid::Grid) = error("implement `get_origin` for $(typeof(grid))")
get_h(grid::Grid) = error("implement `get_h` for $(typeof(grid))")

# ----------------------------
# Derived Utility Methods
# ----------------------------

get_dim(grid::Grid) = length(get_origin(grid))

"""
    get_pos_by_coord(grid::Grid{N, T}, x::SVector{N, T}) -> NTuple{N, Int}

Returns the discrete position (grid indices) corresponding to a coordinate `x`.

- The **cell (0,0) is centered** between `-h/2` and `+h/2`.
- `h` is the length of a grid cell in each dimension.
"""
function get_pos_by_coord(grid::Grid{N}, x) where {N}
    orig = get_origin(grid)
    h = get_h(grid)
    return ntuple(i -> round(Int, (x[i] - orig[i]) / h[i]), Val(N))
end

# The element type is named rather than inferred: a zero-dimensional grid — the input space of a
# hybrid mode whose only control is the switch — has nothing to infer it from, and
# `pos .* get_h(grid)` would yield `SVector{0, Union{}}`. Such a grid holds exactly one cell,
# whose coordinate is the empty vector.
function get_coord_by_pos(grid::Grid{N, T}, pos) where {N, T}
    orig = get_origin(grid)
    h = get_h(grid)
    return SVector{N, T}(ntuple(i -> orig[i] + pos[i] * h[i], Val(N)))
end

# Grid-index bounds of the cells covered by `rect`, one `UnitRange` per
# dimension; an empty range means no cell qualifies in that dimension.
function get_pos_lims_inner(grid::Grid{N}, rect; tol = 1e-6) where {N}
    orig = get_origin(grid)
    h = get_h(grid)
    return ntuple(
        i ->
            ceil(Int, (LazySets.low(rect, i) - tol - orig[i]) / h[i] + 0.5):floor(
                Int,
                (LazySets.high(rect, i) + tol - orig[i]) / h[i] - 0.5,
            ),
        Val(N),
    )
end

function get_pos_lims_outer(grid::Grid{N}, rect; tol = 0.0) where {N}
    orig = get_origin(grid)
    h = get_h(grid)
    return ntuple(
        i ->
            ceil(Int, (LazySets.low(rect, i) + tol - orig[i]) / h[i] - 0.5):floor(
                Int,
                (LazySets.high(rect, i) - tol - orig[i]) / h[i] + 0.5,
            ),
        Val(N),
    )
end

function get_pos_center(grid::Grid{N}, rect; tol = 1e-6) where {N}
    orig = get_origin(grid)
    h = get_h(grid)
    return ntuple(
        i ->
            ceil(Int, (LazySets.low(rect, i) - tol - orig[i]) / h[i]):floor(
                Int,
                (LazySets.high(rect, i) + tol - orig[i]) / h[i],
            ),
        Val(N),
    )
end

function get_pos_lims(grid::Grid, rect, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lims_inner(grid, rect)
    elseif incl_mode == OUTER
        return get_pos_lims_outer(grid, rect)
    else
        return get_pos_center(grid, rect)
    end
end

function get_rec(grid::Grid, pos)
    x = get_coord_by_pos(grid, pos)
    r = get_h(grid) / 2.0
    return LazySets.Hyperrectangle(x, r)
end

get_elem_by_pos(grid::Grid, pos) = get_rec(grid, pos)
get_elem_by_coord(grid::Grid, x) = get_elem_by_pos(grid, get_pos_by_coord(grid, x))
get_all_pos_by_coord(grid::Grid, x) = [get_pos_by_coord(grid, x)]

# Whether cells overlap (a coordinate can belong to several cells). False for a
# tiling grid, true for the ellipsoidal grid whose cells are overlapping balls.
has_overlapping_cells(grid) = false

get_volume(grid::Grid) = prod(get_h(grid))

@recipe function f(grid::Grid, pos; dims = [1, 2])
    @series begin
        UT.project_set(get_elem_by_pos(grid, pos), dims)
    end
end

function get_pos_from_set(grid, rect::LazySets.AbstractHyperrectangle, incl_mode::INCL_MODE)
    return Iterators.product(get_pos_lims(grid, rect, incl_mode)...)
end

# LazySets cannot decide disjointness for every set pair; keeping an
# undecidable candidate is sound for an over-approximation.
function _cell_intersects(cell, S)
    return try
        !LazySets.isdisjoint(cell, S)
    catch
        true
    end
end

# Cells covered by an arbitrary bounded `LazySet`: candidate cells come from
# the bounding-box index ranges, then each candidate is certified per mode —
# OUTER keeps cells intersecting `S` (sound over-approximation), INNER cells
# fully inside (exact for convex `S`, via the cell corners), CENTER cells
# whose center lies in `S`.
function get_pos_from_set(grid, S::LazySets.LazySet, incl_mode::INCL_MODE)
    candidates = Iterators.product(get_pos_lims(grid, UT._outer_box(S), OUTER)...)
    if incl_mode == OUTER
        return (pos for pos in candidates if _cell_intersects(get_rec(grid, pos), S))
    elseif incl_mode == INNER
        return (
            pos for pos in candidates if
            all(v ∈ S for v in LazySets.vertices_list(get_rec(grid, pos)))
        )
    else
        return (pos for pos in candidates if get_coord_by_pos(grid, pos) ∈ S)
    end
end

get_pos_from_set(grid, ::LazySets.EmptySet, incl_mode::INCL_MODE) = ()

function get_pos_from_set(grid, U::LazySets.UnionSetArray, incl_mode::INCL_MODE)
    return Iterators.flatten(get_pos_from_set(grid, s, incl_mode) for s in U.array)
end

function get_pos_from_set(grid::Grid{N}, S::UT.SetMinus, incl_mode::INCL_MODE) where {N}
    inv = invert_incl_mode(incl_mode)
    Bpos = Set{NTuple{N, Int}}()
    for pos in get_pos_from_set(grid, UT.minus_hole(S), inv)
        push!(Bpos, pos)
    end
    return (
        pos for
        pos in get_pos_from_set(grid, UT.minus_included(S), incl_mode) if !(pos in Bpos)
    )
end
# ----------------------------
# Concrete types
# ----------------------------

"""
    GridFree{N,T} <: Grid{N,T}

Uniform grid on unbounded space, centered at `orig` and with steps set by the vector `h`.
"""
struct GridFree{N, T} <: Grid{N, T}
    orig::SVector{N, T}
    h::SVector{N, T}
end

function GridFree(h::SVector{N, T}; zero_origin::Bool = true) where {N, T}
    orig = zero_origin ? zero(h) : h ./ 2
    return GridFree(orig, h)
end

get_origin(grid::GridFree) = grid.orig
get_h(grid::GridFree) = grid.h

"""
    GridEllipsoidalRectangular{N,T} <: Grid{N,T}

Uniform grid on rectangular space `rect`, centered at `orig` and with steps set by the vector `h`.
Cells are (possibly overlapping) ellipsoids defined at each grid point `c` as `(x-c)'P(x-c) ≤ 1`.

`P` (the quadratic form, hot membership loop) and its inverse `Q` (the LazySets
shape matrix, cell construction) are both stored so neither path re-inverts.
"""
struct GridEllipsoidalRectangular{N, T} <: Grid{N, T}
    underlying_grid::GridFree{N, T}
    P::SMatrix{N, N}
    Q::SMatrix{N, N}
end

function GridEllipsoidalRectangular(
    underlying_grid::GridFree{N, T},
    P::SMatrix{N, N},
) where {N, T}
    Q = UT._symmetrize(inv(P))
    return GridEllipsoidalRectangular{N, T}(underlying_grid, P, Q)
end

function GridEllipsoidalRectangular(orig::SVector{N, T}, h::SVector{N, T}, P) where {N, T}
    return GridEllipsoidalRectangular(GridFree(orig, h), SMatrix{N, N}(P))
end

get_origin(grid::GridEllipsoidalRectangular) = get_origin(grid.underlying_grid)
get_h(grid::GridEllipsoidalRectangular) = get_h(grid.underlying_grid)
get_P(grid::GridEllipsoidalRectangular) = grid.P
has_overlapping_cells(grid::GridEllipsoidalRectangular) = true

function get_elem_by_pos(grid::GridEllipsoidalRectangular, pos)
    return LazySets.Ellipsoid(
        collect(get_coord_by_pos(grid, pos)),
        Matrix(grid.Q);
        check_posdef = false,
    )
end

function get_all_pos_by_coord(grid::GridEllipsoidalRectangular{N}, x) where {N}
    center = get_pos_by_coord(grid, x)
    all_pos = typeof(center)[]
    # Neighbouring cells within the {-1,0,1}^N offset box; tuple product avoids
    # the per-call matrix/`eachrow` allocation of `repeat([-1 0 1], N)`.
    for dpos in Iterators.product(ntuple(_ -> (-1, 0, 1), N)...)
        pos = dpos .+ center
        coord = get_coord_by_pos(grid, pos)
        if (x - coord)'grid.P * (x - coord) ≤ 1
            push!(all_pos, pos)
        end
    end
    return all_pos
end
