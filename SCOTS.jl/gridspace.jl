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

function get_pos_lim_inner(grid_space, rect)
    lbI =  ceil.(Int, (rect.lb - grid_space.orig)./grid_space.h .+ 0.5)
    ubI = floor.(Int, (rect.ub - grid_space.orig)./grid_space.h .- 0.5)
    return HyperRectangle(lbI, ubI)
end

function get_pos_lim_outer(grid_space, rect)
    lbI =  ceil.(Int, (rect.lb - grid_space.orig) ./ grid_space.h .- 0.5)
    ubI = floor.(Int, (rect.ub - grid_space.orig) ./ grid_space.h .+ 0.5)
    return HyperRectangle(lbI, ubI)
end

function get_pos_lim(grid_space, rect, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lim_inner(grid_space, rect)
    else
        return get_pos_lim_outer(grid_space, rect)
    end
end

function _make_iterator_from_lims(rect)
    return (collect(pos) for pos in Iterators.product((UnitRange(x...) for x in zip(rect.lb, rect.ub))...))
end

function add_to_gridspace_by_pos!(grid_space::GridSpaceHash, pos)
    push!(grid_space.elems, hash(pos) => pos)
end

function add_to_gridspace_by_coords!(grid_space, x)
    add_to_gridspace_by_pos!(grid_space, get_pos_by_coords(grid_space, x))
end

function add_to_gridspace!(grid_space, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lim(grid_space, rect, incl_mode)
    pos_coll = _make_iterator_from_lims(rectI)
    for pos in pos_coll
        add_to_gridspace_by_pos!(grid_space, pos)
    end
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

function remove_from_gridspace!(grid_space, rect::HyperRectangle, incl_mode::INCL_MODE)
    rectI = get_pos_lim(grid_space, rect, incl_mode)
    pos_coll = _make_iterator_from_lims(rectI)
    if length(pos_coll) < get_gridspace_size(grid_space)
        for pos in pos_coll
            remove_from_gridspace_by_pos!(grid_space, pos)
        end
    else
        for rp in enumerate_gridspace_ref_pos(grid_space)
            if rp[2] in rectI
                remove_from_gridspace_by_ref!(grid_space, rp[1])
            end
        end
    end
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
