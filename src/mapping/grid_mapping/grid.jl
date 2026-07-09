"""
    abstract type Grid{N, T} end

Defines an abstract type for grid-based structures in `N` dimensions with floating-point values `T`.
"""
abstract type Grid{N, T} end

# ----------------------------
# Required API for Grid Domains
# ----------------------------
function get_origin(grid::Grid) end
function get_h(grid::Grid) end

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

function get_coord_by_pos(grid::Grid, pos)
    return get_origin(grid) + pos .* get_h(grid)
end

_ranges(rect::UT.HyperRectangle{N, T}) where {N, T} =
    ntuple(i -> (rect.lb[i]:rect.ub[i]), Val(N))

function get_pos_lims_inner(grid::Grid{N}, rect; tol = 1e-6) where {N}
    orig = get_origin(grid)
    h = get_h(grid)
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - tol - orig[i]) / h[i] + 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] + tol - orig[i]) / h[i] - 0.5), Val(N))
    return UT.HyperRectangle(lbI, ubI)
end

function get_pos_lims_outer(grid::Grid{N}, rect; tol = 0.0) where {N}
    orig = get_origin(grid)
    h = get_h(grid)
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] + tol - orig[i]) / h[i] - 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - tol - orig[i]) / h[i] + 0.5), Val(N))
    return UT.HyperRectangle(lbI, ubI)
end

function get_pos_center(grid::Grid{N}, rect; tol = 1e-6) where {N}
    orig = get_origin(grid)
    h = get_h(grid)
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - tol - orig[i]) / h[i]), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] + tol - orig[i]) / h[i]), Val(N))
    return UT.HyperRectangle(lbI, ubI)
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
    return UT.HyperRectangle(x - r, x + r)
end

get_elem_by_pos(grid::Grid, pos) = get_rec(grid, pos)
get_elem_by_coord(grid::Grid, x) = get_elem_by_pos(grid, get_pos_by_coord(grid, x))
get_all_pos_by_coord(grid::Grid, x) = [get_pos_by_coord(grid, x)]
is_state_cover(grid) = false

function get_volume(grid::Grid)
    r = get_h(grid) / 2.0
    return UT.volume(UT.HyperRectangle(-r, r))
end

@recipe function f(grid::Grid, pos; dims = [1, 2])
    @series begin
        dims := dims
        get_elem_by_pos(grid, pos)
    end
end

function get_pos_from_set(grid, rect::UT.HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(grid, rect, incl_mode)
    return Iterators.product(_ranges(rectI)...)
end

get_pos_from_set(grid, ::UT.EmptyRegion, incl_mode::INCL_MODE) = ()

function get_pos_from_set(grid, U::UT.SetUnion, incl_mode::INCL_MODE)
    return Iterators.flatten(get_pos_from_set(grid, s, incl_mode) for s in U.array)
end

function get_pos_from_set(grid, S::UT.SetMinus, incl_mode::INCL_MODE)
    inv = _invInclMode(incl_mode)
    Bpos = Set{Any}()
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
"""
struct GridEllipsoidalRectangular{N, T} <: Grid{N, T}
    underlying_grid::GridFree{N, T}
    P::SMatrix{N, N}
end

function GridEllipsoidalRectangular(orig::SVector{N, T}, h::SVector{N, T}, P) where {N, T}
    return GridEllipsoidalRectangular{N, T}(GridFree(orig, h), P)
end

get_origin(grid::GridEllipsoidalRectangular) = get_origin(grid.underlying_grid)
get_h(grid::GridEllipsoidalRectangular) = get_h(grid.underlying_grid)
get_P(grid::GridEllipsoidalRectangular) = grid.P
is_state_cover(grid::GridEllipsoidalRectangular) = true

function get_elem_by_pos(grid::GridEllipsoidalRectangular, pos)
    return UT.Ellipsoid(collect(grid.P), collect(get_coord_by_pos(grid, pos)))
end

function get_all_pos_by_coord(grid::GridEllipsoidalRectangular{N}, x) where {N}
    center = get_pos_by_coord(grid, x)
    all_pos = typeof(center)[]
    for dpos in Iterators.product(eachrow(repeat([-1 0 1], N))...)
        coord = get_coord_by_pos(grid, dpos .+ center)
        if (x - coord)'grid.P * (x - coord) ≤ 1
            push!(all_pos, (dpos .+ center))
        end
    end
    return all_pos
end

"""
    DeformedGrid{N, T} <: Grid{N, T}

Represents a deformed version of a `GridFree` grid, where points are mapped via an invertible transformation `f` and its inverse `fi`.

# Fields
- `grid::GridFree{N, T}` : The underlying grid.
- `f::Function` : The forward transformation (physical -> deformed).
- `fi::Function` : The inverse transformation (deformed -> physical).
- `A::Union{Nothing, SMatrix{N, N, T, N*N}}` : Optional linear transformation matrix for volume calculations.
"""
struct DeformedGrid{N, T} <: Grid{N, T}
    underlying_grid::GridFree{N, T}
    f::Function
    fi::Function
    A::Union{Nothing, Any}
end

function DeformedGrid(
    grid::GridFree{N, T},
    f::Function,
    fi::Function;
    A = nothing,
) where {N, T}
    return DeformedGrid{N, T}(grid, f, fi, A)
end

get_origin(grid::DeformedGrid) = grid.f(get_origin(grid.underlying_grid))
get_h(grid::DeformedGrid) = get_h(grid.underlying_grid)

get_pos_by_coord(grid::DeformedGrid, x) = get_pos_by_coord(grid.underlying_grid, grid.fi(x))
get_coord_by_pos(grid::DeformedGrid, pos) =
    grid.f(get_coord_by_pos(grid.underlying_grid, pos))
get_elem_by_pos(grid::DeformedGrid, pos) =
    UT.DeformedRectangle(get_rec(grid.underlying_grid, pos), grid.f)
get_pos_lims(grid::DeformedGrid, rect, incl_mode::INCL_MODE) =
    get_pos_lims(grid.underlying_grid, rect, incl_mode)

function get_volume(grid::DeformedGrid)
    return grid.A !== nothing ? abs(LA.det(grid.A)) * get_volume(grid.underlying_grid) :
           get_volume(grid.underlying_grid)
end
