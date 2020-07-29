function NewTransitionMapHash(X_grid::GridSpaceHash, U_grid::GridSpaceHash, Y_grid::GridSpaceHash)
    elems = Tuple{UInt64,UInt64,UInt64}[]
    infos = Dict("sorted" => true, "unique" => true)
    return TransitionMapHash(X_grid, U_grid, Y_grid, elems, infos)
end

function add_transition_by_ref!(trans_map::TransitionMapHash, x_ref, u_ref, y_ref)
    push!(trans_map.elems, (x_ref, u_ref, y_ref))
    trans_map.infos["sorted"] = false
    trans_map.infos["unique"] = false
end

function _check_sorted_unique(trans_map::TransitionMapHash)
    if ~trans_map.infos["sorted"]
        sort!(trans_map.elems)
        trans_map.infos["sorted"] = true
    end
    if ~trans_map.infos["unique"]
        unique!(trans_map.elems)
        trans_map.infos["unique"] = true
    end
end

function get_transition_image(trans_map::TransitionMapHash, x_ref)
    _check_sorted_unique(trans_map)
    idx_list_x = searchsorted(trans_map.elems, (x_ref,), by = x -> x[1])
    trans_coll_x = trans_map.elems[idx_list_x]
    u_ref_coll = unique(trans[2] for trans in trans_coll_x)
    uy_ref_coll = Vector{Pair{UInt64,Vector{UInt64}}}(undef, length(u_ref_coll))
    for (i, u_ref) in enumerate(u_ref_coll)
        idx_list_u = searchsorted(trans_coll_x, (x_ref, u_ref), by = x -> x[2])
        uy_ref_coll[i] = u_ref => [trans_coll_x[idx][3] for idx in idx_list_u]
    end
    return uy_ref_coll
end

function get_transition_image(trans_map::TransitionMapHash, x_ref, u_ref)
    _check_sorted_unique(trans_map)
    idx_list = searchsorted(trans_map.elems, (x_ref, u_ref), by = x -> x[1:2])
    return [u_ref => [trans_map.elems[idx][3] for idx in idx_list]]
end

function is_xref_controllable(trans_map::TransitionMapHash, x_ref)
    _check_sorted_unique(trans_map)
    idx_list_x = searchsorted(trans_map.elems, Tuple(x_ref), by = x -> x[1])
    return ~isempty(idx_list_x)
end

function remove_all_transitions!(trans_map)
    empty!(trans_map.elems)
    trans_map.infos["sorted"] = true
    trans_map.infos["unique"] = true
end
