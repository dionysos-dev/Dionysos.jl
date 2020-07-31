function NewSubSet(grid_space::GridSpaceHash{N}) where N
    return SubSetHash{N}(grid_space, UInt64[], true, true)
end

function ensure_sorted!(sub_set::SubSetHash)
    if ~sub_set.issorted
        # display("subset not sorted")
        sort!(sub_set.elems)
        sub_set.issorted = true
    end
end

function ensure_unique!(sub_set::SubSetHash)
    if ~sub_set.isunique
        # display("subset not unique")
        unique!(sub_set.elems)
        sub_set.isunique = true
    end
end

function sizehint_subset!(sub_set::SubSetHash, size_max)
    sizehint!(sub_set.elems, size_max)
end

function add_to_subset_by_ref!(sub_set::SubSetHash, ref)
    push!(sub_set.elems, ref)
    sub_set.issorted = false
    sub_set.isunique = false
end

function add_to_subset_by_new_ref!(sub_set::SubSetHash, ref)
    push!(sub_set.elems, ref)
    sub_set.issorted = false
end

function add_to_subset_by_ref_coll!(sub_set::SubSetHash, ref_coll)
    append!(sub_set.elems, ref_coll)
    if ~isempty(ref_coll)
        sub_set.issorted = false
        sub_set.isunique = false
    end
end

function add_to_subset_by_new_ref_coll!(sub_set::SubSetHash, ref_coll)
    append!(sub_set.elems, ref_coll)
    if ~isempty(ref_coll)
        sub_set.issorted = false
    end
end

function union_subsets!(sub_set1::SubSetHash, sub_set2::SubSetHash)
    append!(sub_set1.elems, sub_set2.elems)
    if ~is_subset_empty(sub_set2)
        sub_set1.issorted = false
        sub_set1.isunique = false
    end
end

function add_to_subset_all!(sub_set::SubSetHash)
    add_to_subset_by_ref_coll!(sub_set, enumerate_gridspace_ref(sub_set.grid_space))
end

function add_to_subset_by_pos!(sub_set, pos)
    ref = get_ref_by_pos(sub_set.grid_space, pos)
    if ref !== sub_set.grid_space.overflow_ref
        add_to_subset_by_ref!(sub_set, ref)
    end
end

function add_to_subset_by_pos_coll!(sub_set, pos_coll)
    for pos in pos_coll
        add_to_subset_by_pos!(sub_set, pos)
    end
end

function add_to_subset_by_coords!(sub_set, x)
    add_to_subset_by_pos!(sub_set, get_pos_by_coords(sub_set.grid_space, x))
end

function add_to_subset!(sub_set, rect, incl_mode::INCL_MODE)
    rectI = get_pos_lims(sub_set.grid_space, rect, incl_mode)
    pos_iter = _make_iterator_from_lims(rectI)
    if length(pos_iter) < get_gridspace_size(sub_set.grid_space)
        add_to_subset_by_pos_coll!(sub_set, pos_iter)
    else
        for rp in enumerate_gridspace_ref_pos(sub_set.grid_space)
            if rp[2] in rectI
                add_to_subset_by_ref!(sub_set, rp[1])
            end
        end
    end
end

function remove_from_subset_by_ref!(sub_set::SubSetHash, ref)
    setdiff!(sub_set.elems, ref)
end

function remove_from_subset_by_ref_coll!(sub_set::SubSetHash, ref_coll)
    setdiff!(sub_set.elems, ref_coll)
end

function setdiff_subsets!(sub_set1::SubSetHash, sub_set2::SubSetHash)
    setdiff!(sub_set1.elems, sub_set2.elems)
end

function remove_from_subset_all!(sub_set::SubSetHash)
    empty!(sub_set.elems)
    sub_set.issorted = true
    sub_set.isunique = true
end

function remove_from_subset_by_pos!(sub_set, pos)
    ref = get_ref_by_pos(sub_set.grid_space, pos)
    if ref !== sub_set.grid_space.overflow_ref
        remove_from_subset_by_ref!(sub_set, ref)
    end
end

function remove_from_subset_by_pos_coll!(sub_set, pos_coll)
    for pos in pos_coll
        remove_from_subset_by_pos!(sub_set, pos)
    end
end

function remove_from_subset_by_coords!(sub_set, x)
    remove_from_subset_by_pos!(sub_set, get_pos_by_coords(sub_set.grid_space, x))
end

function remove_from_subset!(sub_set, rect, incl_mode::INCL_MODE)
    rectI = get_pos_lims(sub_set.grid_space, rect, incl_mode)
    pos_iter = _make_iterator_from_lims(rectI)
    if length(pos_iter) < get_subset_size(sub_set)
        remove_from_subset_by_pos_coll!(sub_set, pos_iter)
    else
        for ref in enumerate_subset_ref(sub_set)
            if ref !== sub_set.grid_space.overflow_ref
                pos = get_pos_by_ref(sub_set.grid_space, ref)
                if pos in rectI
                    remove_from_subset_by_ref!(sub_set, ref)
                end
            end
        end
    end
end

function enumerate_subset_ref(sub_set::SubSetHash)
    ensure_unique!(sub_set)
    return sub_set.elems
end

function is_ref_in_subset(sub_set::SubSetHash, ref)
    ensure_sorted!(sub_set)
    ensure_unique!(sub_set)
    return ~isempty(searchsorted(sub_set.elems, ref))
end

function get_subset_size(sub_set::SubSetHash)
    ensure_unique!(sub_set)
    return length(sub_set.elems)
end

function is_subset_empty(sub_set::SubSetHash)
    return isempty(sub_set.elems)
end

function sizehint_subset!(sub_set::SubSetHash, size_max)
    sizehint!(sub_set.elems, size_max)
end
