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
    uy_ref_coll = Pair{UInt64,Vector{UInt64}}[]
    if isempty(idx_list_x)
        return uy_ref_coll
    end
    cur_i = idx_list_x[1]
    cur_u = uy_ref_coll[i][1]
    for i in idx_list_x[2:end]
        _, u, y = trans_map.elems[i]
        if cur_u != u
            push!(uy_ref_coll, cur_u => UInt64[trans_map.elems[j][3] for j in cur_i:(i-1)])
            cur_u = u
            cur_i = i
        end
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
