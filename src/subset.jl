abstract type SubSet{N} end

struct SubSetList{N} <: SubSet{N}
    gridspace::GridSpace{N}
    elems::Set{NTuple{N, Int}}
end

function NewSubSet(gridspace::GridSpaceList{N}) where N
    return SubSetList{N}(gridspace, Set{NTuple{N, Int}}())
end

function add_pos!(subset::SubSetList, pos)
    if pos ∈ subset.gridspace
        push!(subset.elems, pos)
    end
end

function Base.union!(subset1::SubSetList, subset2::SubSetList)
    union!(subset1.elems, subset2.elems)
end

function add_all!(subset::SubSet)
    for pos in enum_pos(subset.gridspace)
        add_pos!(subset, pos)
    end
end

function add_coord!(subset::SubSet, x)
    add_pos!(subset, get_pos_by_coord(subset.gridspace, x))
end

function add_set!(subset::SubSet, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(subset.gridspace, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(subset.gridspace)
        for pos in pos_iter
            add_pos!(subset, pos)
        end
    else
        for pos in Iterators.filter(x -> x in rectI, enum_pos(subset.gridspace))
            add_pos!(subset, pos)
        end
    end
end

function remove_pos!(subset::SubSetList, pos)
    delete!(subset.elems, pos)
end

function Base.setdiff!(subset1::SubSetList, subset2::SubSetList)
    setdiff!(subset1.elems, subset2.elems)
end

function Base.empty!(subset::SubSetList)
    empty!(subset.elems)
end

function remove_coord!(subset::SubSet, x)
    remove_pos!(subset, get_pos_by_coord(subset.gridspace, x))
end

function remove_set!(subset::SubSet, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lims(subset.gridspace, rect, incl_mode)
    pos_iter = Iterators.product(_ranges(rectI)...)
    if length(pos_iter) < get_ncells(subset)
        for pos in pos_iter
            remove_pos!(subset, pos)
        end
    else
        for pos in Iterators.filter(x -> x in rectI, enum_pos(subset.gridspace))
            remove_pos!(subset, pos)
        end
    end
end

function enum_pos(subset::SubSetList)
    return subset.elems
end

function Base.in(pos, subset::SubSetList)
    return in(pos, subset.elems)
end

function get_ncells(subset::SubSetList)
    return length(subset.elems)
end

function Base.isempty(subset::SubSetList)
    return isempty(subset.elems)
end

function Base.issubset(subset1::SubSetList, subset2::SubSetList)
    return issubset(subset1.elems, subset2.elems)
end

function get_somepos(subset::SubSetList)
    return first(subset.elems)
end
