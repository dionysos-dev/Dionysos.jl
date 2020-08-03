abstract type GridSpace{N} end

struct GridSpaceList{N} <: GridSpace{N}
    orig::NTuple{N, Float64}
    h::NTuple{N, Float64}
    elems::Set{NTuple{N, Int}}
end

function NewGridSpaceList(orig::NTuple{N, Float64}, h::NTuple{N, Float64}) where N
    return GridSpaceList{N}(orig, h, D, Set{NTuple{N, Int}}())
end

function get_pos_by_coord(gridspace, x::Tuple)
    return round.(Int, (x .- gridspace.orig)./gridspace.h)
end

function get_coord_by_pos(gridspace, pos)
    return gridspace.orig .+ pos.*gridspace.h
end

function get_pos_lims_inner(gridspace, rect)
    lbI = ceil.(Int, (rect.lb .- gridspace.orig)./gridspace.h .+ 0.5)
    ubI = floor.(Int, (rect.ub .- gridspace.orig)./gridspace.h .- 0.5)
    return HyperRectangle(lbI, ubI)
end

function get_pos_lims_outer(gridspace, rect)
    lbI = ceil.(Int, (rect.lb .- gridspace.orig)./gridspace.h .- 0.5)
    ubI = floor.(Int, (rect.ub .- gridspace.orig)./gridspace.h .+ 0.5)
    return HyperRectangle(lbI, ubI)
end

function get_pos_lims(gridspace, rect, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lims_inner(gridspace, rect)
    else
        return get_pos_lims_outer(gridspace, rect)
    end
end

function _ranges(rect::HyperRectangle{NTuple{N, T}}) where {T, N}
    return ntuple(i -> UnitRange(rect.lb[i], rect.ub[i]), Val(N))
end

function add_position!(gridspace::GridSpaceList, pos)
    push!(gridspace.elems, pos)
end

function add_coord!(gridspace, x)
    add_pos!(gridspace, get_pos_by_coord(gridspace, x))
end

function add_set!(gridspace, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(gridspace, rect, incl_mode)
    for pos in Iterators.product(_ranges(rectI)...)
        add_pos!(gridspace, pos)
    end
end

function remove_pos!(gridspace::GridSpaceList, pos)
    delete!(gridspace.elems, pos)
end

function remove_coord!(gridspace, x)
    remove_pos!(gridspace, get_pos_by_coord(gridspace, x))
end

function remove_set!(gridspace, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(gridspace, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_gridspace_size(gridspace)
        for pos in pos_iter
            remove_pos!(gridspace, pos)
        end
    else
        for pos in Iterators.filter(x -> x in rectI, enum_pos(gridspace))
            remove_pos!(gridspace, pos)
        end
    end
end

function has_pos(gridspace::GridSpaceList, pos)
    return in(pos, gridspace.elems)
end

function enum_pos(gridspace::GridSpaceList)
    return gridspace.elems
end

function get_ncells(gridspace::GridSpaceList)
    return length(gridspace.elems)
end

function isempty(gridspace::GridSpaceList)
    return isempty(gridspace.elems)
end
