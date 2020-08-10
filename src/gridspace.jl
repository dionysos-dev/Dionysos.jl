abstract type GridSpace{N,T} end

struct GridSpaceList{N,T} <: GridSpace{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
    elems::Set{NTuple{N,Int}}
end

function NewGridSpaceList(orig::SVector{N,T}, h::SVector{N,T}) where {N,T}
    return GridSpaceList(orig, h, Set{NTuple{N,Int}}())
end

function get_pos_by_coord(gridspace::GridSpace{N}, x) where N
    return ntuple(i -> round(Int, (x[i] - gridspace.orig[i])/gridspace.h[i]), Val(N))
end

function get_coord_by_pos(gridspace, pos)
    return gridspace.orig .+ pos.*gridspace.h
end

function get_pos_lims_inner(gridspace::GridSpace{N}, rect) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - gridspace.orig[i])/gridspace.h[i] + 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - gridspace.orig[i])/gridspace.h[i] - 0.5), Val(N))
    return HyperRectangle(lbI, ubI)
end

function get_pos_lims_outer(gridspace::GridSpace{N}, rect) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - gridspace.orig[i])/gridspace.h[i] - 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - gridspace.orig[i])/gridspace.h[i] + 0.5), Val(N))
    return HyperRectangle(lbI, ubI)
end

function get_pos_lims(gridspace, rect, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lims_inner(gridspace, rect)
    else
        return get_pos_lims_outer(gridspace, rect)
    end
end

function _ranges(rect::HyperRectangle{NTuple{N,T}}) where {N,T}
    return ntuple(i -> UnitRange(rect.lb[i], rect.ub[i]), Val(N))
end

function add_pos!(gridspace::GridSpaceList, pos)
    push!(gridspace.elems, pos)
end

function add_coord!(gridspace::GridSpace, x)
    add_pos!(gridspace, get_pos_by_coord(gridspace, x))
end

function add_set!(gridspace::GridSpace, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(gridspace, rect, incl_mode)
    for pos in Iterators.product(_ranges(rectI)...)
        add_pos!(gridspace, pos)
    end
end

function remove_pos!(gridspace::GridSpaceList, pos)
    delete!(gridspace.elems, pos)
end

function remove_coord!(gridspace::GridSpace, x)
    remove_pos!(gridspace, get_pos_by_coord(gridspace, x))
end

function remove_set!(gridspace::GridSpace, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(gridspace, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(gridspace)
        for pos in pos_iter
            remove_pos!(gridspace, pos)
        end
    else
        for pos in enum_pos(gridspace)
            if pos ∈ rectI
                remove_pos!(gridspace, pos)
            end
        end
    end
end

function Base.in(pos, gridspace::GridSpaceList)
    return in(pos, gridspace.elems)
end

function enum_pos(gridspace::GridSpaceList)
    return gridspace.elems
end

function get_ncells(gridspace::GridSpaceList)
    return length(gridspace.elems)
end

function Base.isempty(gridspace::GridSpaceList)
    return isempty(gridspace.elems)
end

function get_somepos(gridspace::GridSpaceList)
    return first(gridspace.elems)
end
