function NewGridSpaceHash(orig, h)
    @assert length(orig) == length(h)
    D = Dict{UInt64,Vector{Int}}()
    overflow_pos = Vector{Int}(undef, length(orig))
    return GridSpaceHash(length(orig), orig, h, D, UInt64(0), overflow_pos)
end

function get_pos_by_coords(grid_space, x)
    return round.(Int, (x - grid_space.orig)./grid_space.h)
end

function get_coords_by_pos(grid_space, pos)
    return grid_space.orig + pos.*grid_space.h
end

function get_pos_lim_from_box_inner(grid_space, lb, ub)
    lbI = ceil.(Int, (lb - grid_space.orig)./grid_space.h .+ 0.5)
    ubI = floor.(Int, (ub - grid_space.orig)./grid_space.h .- 0.5)
    return (lbI, ubI)
end

function get_pos_lim_from_box_outer(grid_space, lb, ub)
    lbI = ceil.(Int, (lb - grid_space.orig)./grid_space.h .- 0.5)
    ubI = floor.(Int, (ub - grid_space.orig)./grid_space.h .+ 0.5)
    return (lbI, ubI)
end

function get_pos_lim_from_box(grid_space, lb, ub, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lim_from_box_inner(grid_space, lb, ub)
    else
        return get_pos_lim_from_box_outer(grid_space, lb, ub)
    end
end

function _make_iterator_from_lims(lbI, ubI)
    return (collect(pos) for pos in Iterators.product((UnitRange(x...) for x in zip(lbI, ubI))...))
end

function _is_pos_in_lims(pos, lbI, ubI)
    return all(lbI .<= pos .<= ubI)
end

function add_to_gridspace_by_pos!(grid_space::GridSpaceHash, pos)
    push!(grid_space.elems, hash(pos) => pos)
end

function add_to_gridspace_by_coords!(grid_space, x)
    add_to_gridspace_by_pos!(grid_space, get_pos_by_coords(grid_space, x))
end

function add_to_gridspace_by_box!(grid_space, lb, ub, incl_mode::INCL_MODE)
    lbI, ubI = get_pos_lim_from_box(grid_space, lb, ub, incl_mode)
    pos_coll = _make_iterator_from_lims(lbI, ubI)
    for pos in pos_coll
        add_to_gridspace_by_pos!(grid_space, pos)
    end
end
function add_to_gridspace!(grid_space, rect::HyperRectangle, incl_mode::INCL_MODE)
    add_to_gridspace_by_box!(grid_space, rect.lb, rect.ub, incl_mode)
end

function remove_from_gridspace_by_ref!(grid_space::GridSpaceHash, ref)
    delete!(grid_space.elems, ref)
end

function remove_from_gridspace_by_pos!(grid_space::GridSpaceHash, pos)
    delete!(grid_space.elems, hash(pos))
end

function remove_from_gridspace_by_coords!(grid_space, x)
    remove_from_gridspace_by_pos!(grid_space, get_pos_by_coords(grid_space, x))
end

function remove_from_gridspace_by_box!(grid_space, lb, ub, incl_mode::INCL_MODE)
    lbI, ubI = get_pos_lim_from_box(grid_space, lb, ub, incl_mode)
    pos_coll = _make_iterator_from_lims(lbI, ubI)
    if length(pos_coll) < get_gridspace_size(grid_space)
        for pos in pos_coll
            remove_from_gridspace_by_pos!(grid_space, pos)
        end
    else
        for rp in enumerate_gridspace_ref_pos(grid_space)
            if _is_pos_in_lims(rp[2], lbI, ubI)
                remove_from_gridspace_by_ref!(grid_space, rp[1])
            end
        end
    end
end

function remove_from_gridspace!(grid_space, rect::HyperRectangle, incl_mode::INCL_MODE)
    remove_from_gridspace_by_box!(grid_space, rect.lb, rect.ub, incl_mode)
end

function get_ref_by_pos(grid_space::GridSpaceHash, pos)
    return getkey(grid_space.elems, hash(pos), grid_space.overflow_ref)
end

function get_pos_by_ref(grid_space::GridSpaceHash, ref)
    return get(grid_space.elems, ref, grid_space.overflow_pos)
end

function get_coords_by_ref(grid_space, ref)
    return get_coords_by_pos(grid_space, get_pos_by_ref(grid_space, ref))
end

function get_ref_by_coords(grid_space, x)
    return get_ref_by_pos(grid_space, get_pos_by_coords(grid_space, x))
end

function enumerate_gridspace_ref(grid_space::GridSpaceHash)
    return keys(grid_space.elems)
end

function enumerate_gridspace_pos(grid_space::GridSpaceHash)
    return values(grid_space.elems)
end

function enumerate_gridspace_ref_pos(grid_space::GridSpaceHash)
    return pairs(grid_space.elems)
end

function get_gridspace_size(grid_space::GridSpaceHash)
    return length(grid_space.elems)
end

function is_gridspace_empty(grid_space::GridSpaceHash)
    return isempty(grid_space.elems)
end
