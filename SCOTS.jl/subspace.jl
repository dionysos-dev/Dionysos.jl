function NewSubSpace(grid_space::GridSpaceHash)
    return SubSpaceHash(grid_space, Set{UInt64}())
end

function add_to_subspace_by_ref!(sub_space::SubSpaceHash, ref)
    push!(sub_space.elems, ref)
end

function add_to_subspace_all!(sub_space::SubSpaceHash)
    for ref in enumerate_gridspace_ref(sub_space.grid_space)
        add_to_subspace_by_ref!(sub_space, ref)
    end
end

function add_to_subspace_by_pos!(sub_space, pos)
    ref = get_ref_by_pos(sub_space.grid_space, pos)
    if ref !== sub_space.grid_space.overflow_ref
        add_to_subspace_by_ref!(sub_space, ref)
    end
end

function add_to_subspace_by_coords!(sub_space, x)
    add_to_subspace_by_pos!(sub_space, get_pos_by_coords(sub_space.grid_space, x))
end

function add_to_subspace_by_box!(sub_space, lb, ub, incl_mode::INCL_MODE)
    lbI, ubI = get_pos_lim_from_box(sub_space.grid_space, lb, ub, incl_mode)
    pos_coll = _make_iterator_from_lims(lbI, ubI)
    if length(pos_coll) < get_gridspace_size(sub_space.grid_space)
        for pos in pos_coll
            add_to_subspace_by_pos!(sub_space, pos)
        end
    else
        for rp in enumerate_gridspace_ref_pos(sub_space.grid_space)
            if _is_pos_in_lims(rp[2], lbI, ubI)
                add_to_subspace_by_ref!(sub_space, rp[1])
            end
        end
    end
end

function add_to_subspace!(sub_space, rect::HyperRectangle, incl_mode::INCL_MODE)
    add_to_subspace_by_box!(sub_space, rect.lb, rect.ub, incl_mode)
end

function remove_from_subspace_by_ref!(sub_space::SubSpaceHash, ref)
    delete!(sub_space.elems, ref)
end

function remove_from_subspace_by_pos!(sub_space, pos)
    delete!(sub_space.elems, get_ref_by_pos(sub_space.grid_space, pos))
end

function remove_from_subspace_by_coords!(sub_space, x)
    remove_from_subspace_by_pos!(sub_space, get_pos_by_coords(sub_space.grid_space, x))
end

function remove_from_subspace_by_box!(sub_space, lb, ub, incl_mode::INCL_MODE)
    lbI, ubI = get_pos_lim_from_box(sub_space.grid_space, lb, ub, incl_mode)
    pos_coll = _make_iterator_from_lims(lbI, ubI)
    if length(pos_coll) < get_subspace_size(sub_space)
        for pos in pos_coll
            remove_from_subspace_by_pos!(sub_space, pos)
        end
    else
        for ref in enumerate_subspace_ref(sub_space)
            if ref !== sub_space.grid_space.overflow_ref
                pos = get_pos_by_ref(sub_space.grid_space, ref)
                if _is_pos_in_lims(pos, lbI, ubI)
                    remove_from_subspace_by_ref!(sub_space, ref)
                end
            end
        end
    end
end

function enumerate_subspace_ref(sub_space::SubSpaceHash)
    return sub_space.elems
end

function is_ref_in_subspace(sub_space::SubSpaceHash, ref)
    return ∈(ref, sub_space.elems)
end

function get_subspace_size(sub_space::SubSpaceHash)
    return length(sub_space.elems)
end

function is_subspace_empty(sub_space::SubSpaceHash)
    return isempty(sub_space.elems)
end
