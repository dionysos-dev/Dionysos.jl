function NewSubSpace(grid_space::GridSpaceHash{N}) where N
    return SubSpaceHash{N}(grid_space, UInt64[], true, true)
end

function ensure_sorted!(sub_space::SubSpaceHash)
    if ~sub_space.issorted
        # display("subspace not sorted")
        sort!(sub_space.elems)
        sub_space.issorted = true
    end
end

function ensure_unique!(sub_space::SubSpaceHash)
    if ~sub_space.isunique
        # display("subspace not unique")
        unique!(sub_space.elems)
        sub_space.isunique = true
    end
end

function sizehint_subspace!(sub_space::SubSpaceHash, size_max)
    sizehint!(sub_space.elems, size_max)
end

function add_to_subspace_by_ref!(sub_space::SubSpaceHash, ref)
    push!(sub_space.elems, ref)
    sub_space.issorted = false
    sub_space.isunique = false
end

function add_to_subspace_by_new_ref!(sub_space::SubSpaceHash, ref)
    push!(sub_space.elems, ref)
    sub_space.issorted = false
end

function add_to_subspace_by_ref_coll!(sub_space::SubSpaceHash, ref_coll)
    append!(sub_space.elems, ref_coll)
    if ~isempty(ref_coll)
        sub_space.issorted = false
        sub_space.isunique = false
    end
end

function add_to_subspace_by_new_ref_coll!(sub_space::SubSpaceHash, ref_coll)
    append!(sub_space.elems, ref_coll)
    if ~isempty(ref_coll)
        sub_space.issorted = false
    end
end

function union_subspaces!(sub_space1::SubSpaceHash, sub_space2::SubSpaceHash)
    append!(sub_space1.elems, sub_space2.elems)
    if ~is_subspace_empty(sub_space2)
        sub_space1.issorted = false
        sub_space1.isunique = false
    end
end

function add_to_subspace_all!(sub_space::SubSpaceHash)
    add_to_subspace_by_ref_coll!(sub_space, enumerate_gridspace_ref(sub_space.grid_space))
end

function add_to_subspace_by_pos!(sub_space, pos)
    ref = get_ref_by_pos(sub_space.grid_space, pos)
    if ref !== sub_space.grid_space.overflow_ref
        add_to_subspace_by_ref!(sub_space, ref)
    end
end

function add_to_subspace_by_pos_coll!(sub_space, pos_coll)
    for pos in pos_coll
        add_to_subspace_by_pos!(sub_space, pos)
    end
end

function add_to_subspace_by_coords!(sub_space, x)
    add_to_subspace_by_pos!(sub_space, get_pos_by_coords(sub_space.grid_space, x))
end

function add_to_subspace_by_box!(sub_space, lb, ub, incl_mode::INCL_MODE)
    lbI, ubI = get_pos_lims_from_box(sub_space.grid_space, lb, ub, incl_mode)
    pos_iter = _make_iterator_from_lims(lbI, ubI)
    if length(pos_iter) < get_gridspace_size(sub_space.grid_space)
        add_to_subspace_by_pos_coll!(sub_space, pos_iter)
    else
        for rp in enumerate_gridspace_ref_pos(sub_space.grid_space)
            if is_pos_in_lims(rp[2], lbI, ubI)
                add_to_subspace_by_ref!(sub_space, rp[1])
            end
        end
    end
end

function remove_from_subspace_by_ref!(sub_space::SubSpaceHash, ref)
    setdiff!(sub_space.elems, ref)
end

function remove_from_subspace_by_ref_coll!(sub_space::SubSpaceHash, ref_coll)
    setdiff!(sub_space.elems, ref_coll)
end

function setdiff_subspaces!(sub_space1::SubSpaceHash, sub_space2::SubSpaceHash)
    setdiff!(sub_space1.elems, sub_space2.elems)
end

function remove_from_subspace_all!(sub_space::SubSpaceHash)
    empty!(sub_space.elems)
    sub_space.issorted = true
    sub_space.isunique = true
end

function remove_from_subspace_by_pos!(sub_space, pos)
    ref = get_ref_by_pos(sub_space.grid_space, pos)
    if ref !== sub_space.grid_space.overflow_ref
        remove_from_subspace_by_ref!(sub_space, ref)
    end
end

function remove_from_subspace_by_pos_coll!(sub_space, pos_coll)
    for pos in pos_coll
        remove_from_subspace_by_pos!(sub_space, pos)
    end
end

function remove_from_subspace_by_coords!(sub_space, x)
    remove_from_subspace_by_pos!(sub_space, get_pos_by_coords(sub_space.grid_space, x))
end

function remove_from_subspace_by_box!(sub_space, lb, ub, incl_mode::INCL_MODE)
    lbI, ubI = get_pos_lims_from_box(sub_space.grid_space, lb, ub, incl_mode)
    pos_iter = _make_iterator_from_lims(lbI, ubI)
    if length(pos_iter) < get_subspace_size(sub_space)
        remove_from_subspace_by_pos_coll!(sub_space, pos_iter)
    else
        for ref in enumerate_subspace_ref(sub_space)
            if ref !== sub_space.grid_space.overflow_ref
                pos = get_pos_by_ref(sub_space.grid_space, ref)
                if is_pos_in_lims(pos, lbI, ubI)
                    remove_from_subspace_by_ref!(sub_space, ref)
                end
            end
        end
    end
end

function enumerate_subspace_ref(sub_space::SubSpaceHash)
    ensure_unique!(sub_space)
    return sub_space.elems
end

function is_ref_in_subspace(sub_space::SubSpaceHash, ref)
    ensure_sorted!(sub_space)
    ensure_unique!(sub_space)
    return ~isempty(searchsorted(sub_space.elems, ref))
end

function get_subspace_size(sub_space::SubSpaceHash)
    ensure_unique!(sub_space)
    return length(sub_space.elems)
end

function is_subspace_empty(sub_space::SubSpaceHash)
    return isempty(sub_space.elems)
end

function sizehint_subspace!(sub_space::SubSpaceHash, size_max)
    sizehint!(sub_space.elems, size_max)
end
