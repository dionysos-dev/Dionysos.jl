abstract type Grid{N,T} end


# Free because later: maybe put bounds lb and ub (e.g., for BDDs)
struct GridFree{N,T} <: Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
end

struct GridRectangular{N,T} <: Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
    rect::Any
end


struct GridEllipsoidalRectangular{N,T} <: Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
    P::SMatrix{N,N}
    rect::Any
end


function get_pos_by_coord(grid::Grid{N}, x) where N
    return ntuple(i -> round(Int, (x[i] - grid.orig[i])/grid.h[i]), Val(N))
end

function get_all_pos_by_coord(grid::GridEllipsoidalRectangular{N}, x) where N
    center = get_pos_by_coord(grid,x)
    all_pos = typeof(center)[]
    for dpos in Iterators.product(eachrow(repeat([-1 0 1],N))...)
        coord = get_coord_by_pos(grid, dpos.+center)
        if (x-coord)'grid.P*(x-coord) ≤ 1
            push!(all_pos, (dpos.+center))
        end
    end
    return all_pos
end


function get_coord_by_pos(grid, pos)
    return grid.orig + pos.*grid.h
end

function get_pos_lims_inner(grid::Grid{N}, rect; tol=10^(-6)) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - tol - grid.orig[i])/grid.h[i] + 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] + tol - grid.orig[i])/grid.h[i] - 0.5), Val(N))
    return HyperRectangle(lbI, ubI)
end

function get_pos_lims_outer(grid::Grid{N}, rect) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - grid.orig[i])/grid.h[i] - 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - grid.orig[i])/grid.h[i] + 0.5), Val(N))
    return UT.HyperRectangle(lbI, ubI)
end

function get_pos_lims(grid, rect, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lims_inner(grid, rect)
    else
        return get_pos_lims_outer(grid, rect)
    end
end

function _ranges(rect::UT.HyperRectangle{NTuple{N,T}}) where {N,T}
    return ntuple(i -> UnitRange(rect.lb[i], rect.ub[i]), Val(N))
end
