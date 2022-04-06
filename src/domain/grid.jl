abstract type Grid{N,T} end


# Free because later: maybe put bounds lb and ub (e.g., for BDDs)
struct GridFree{N,T} <: Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
end

function get_pos_by_coord(grid::Grid{N}, x) where N
    return ntuple(i -> round(Int, (x[i] - grid.orig[i])/grid.h[i]), Val(N))
end

function get_coord_by_pos(grid, pos)
    return grid.orig + pos.*grid.h
end

function get_pos_lims_inner(grid::Grid{N}, rect) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - grid.orig[i])/grid.h[i] + 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - grid.orig[i])/grid.h[i] - 0.5), Val(N))
    return UT.HyperRectangle(lbI, ubI)
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
