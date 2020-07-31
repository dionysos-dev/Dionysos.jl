function NewSymbolicModelHash(X_grid::GridSpaceHash, U_grid::GridSpaceHash, Y_grid::GridSpaceHash)
    elems = Tuple{UInt64,UInt64,UInt64}[]
    return SymbolicModelHash(X_grid, U_grid, Y_grid, elems, true, true)
end

function ensure_sorted!(sym_model::SymbolicModelHash)
    if !sym_model.issorted
        # display("sym_model not sorted")
        sort!(sym_model.elems)
        sym_model.issorted = true
    end
end

function ensure_unique!(sym_model::SymbolicModelHash)
    if !sym_model.isunique
        # display("sym_model not unique")
        unique!(sym_model.elems)
        sym_model.isunique = true
    end
end

function add_to_symmodel_by_refs!(sym_model::SymbolicModelHash, x_ref, u_ref, y_ref)
    push!(sym_model.elems, (x_ref, u_ref, y_ref))
    sym_model.issorted = false
    sym_model.isunique = false
end

function add_to_symmodel_by_new_refs!(sym_model::SymbolicModelHash, x_ref, u_ref, y_ref)
    push!(sym_model.elems, (x_ref, u_ref, y_ref))
    sym_model.issorted = false
end

function add_to_symmodel_by_refs_coll!(sym_model::SymbolicModelHash, refs_coll)
    append!(sym_model.elems, refs_coll)
    if !isempty(refs_coll)
        sym_model.issorted = false
        sym_model.isunique = false
    end
end

function add_to_symmodel_by_new_refs_coll!(sym_model::SymbolicModelHash, refs_coll)
    append!(sym_model.elems, refs_coll)
    if !isempty(refs_coll)
        sym_model.issorted = false
    end
end

function remove_from_symmodel_all!(sym_model)
    empty!(sym_model.elems)
    sym_model.issorted = true
    sym_model.isunique = true
end

# function add_images_by_xref_uref!(Y_sub, sym_model::SymbolicModelHash, x_ref, u_ref)
#     ensure_sorted!(sym_model)
#     ensure_unique!(sym_model)
#     idx_iter = searchsorted(sym_model.elems, (x_ref, u_ref), by = x -> x[1:2])
#     add_to_subset_by_ref_coll!(Y_sub, sym_model.elems[idx][3] for idx in idx_iter)
# end

function add_images_by_xref_uref!(yref_coll, sym_model::SymbolicModelHash, x_ref, u_ref)
    ensure_sorted!(sym_model)
    ensure_unique!(sym_model)
    idx_iter = searchsorted(sym_model.elems, (x_ref, u_ref), by = x -> x[1:2])
    for idx in idx_iter
        push!(yref_coll, sym_model.elems[idx][3])
    end
end

# "Set" (vs "add") assumes Y_sub is empty initially... May fail if not respected
# function set_images_by_xref_uref!(Y_sub::SubSetHash, sym_model::SymbolicModelHash, x_ref, u_ref)
#     add_images_by_xref_uref!(Y_sub, sym_model, x_ref, u_ref)
#     Y_sub.issorted = true
#     Y_sub.isunique = true
# end

function is_xref_controllable(sym_model::SymbolicModelHash, x_ref)
    ensure_sorted!(sym_model)
    idx_iter_x = searchsorted(sym_model.elems, (x_ref,), by = x -> x[1])
    return !isempty(idx_iter_x)
end

# function add_inputs_by_xref_ysub!(U_sub, sym_model::SymbolicModelHash, x_ref, Y_sub)
#     ensure_sorted!(sym_model)
#     ensure_unique!(sym_model)
#     idx_iter = searchsorted(sym_model.elems, (x_ref,), by = x -> x[1])
#     uref_prev = U_sub.grid_space.overflow_ref
#     refs = (x_ref, uref_prev, Y_sub.grid_space.overflow_ref)
#     all_in = false
#
#     for idx in idx_iter
#         refs = sym_model.elems[idx]
#         if refs[2] != uref_prev
#             if all_in
#                 add_to_subset_by_ref!(U_sub, uref_prev)
#             else
#                 all_in = true
#             end
#             uref_prev = refs[2]
#         end
#         all_in = all_in && is_ref_in_subset(Y_sub, refs[3])
#     end
#
#     if all_in
#         add_to_subset_by_ref!(U_sub, refs[2])
#     end
# end

function add_inputs_by_xref_ysub!(uref_coll, sym_model::SymbolicModelHash, x_ref, Y_sub)
    ensure_sorted!(sym_model)
    ensure_unique!(sym_model)
    idx_iter = searchsorted(sym_model.elems, (x_ref,), by = x -> x[1])
    uref_prev = sym_model.U_grid.overflow_ref
    refs = (x_ref, uref_prev, sym_model.Y_grid.overflow_ref)
    all_in = false

    for idx in idx_iter
        refs = sym_model.elems[idx]
        if refs[2] != uref_prev
            if all_in
                push!(uref_coll, uref_prev)
            else
                all_in = true
            end
            uref_prev = refs[2]
        end
        all_in = all_in && is_ref_in_subset(Y_sub, refs[3])
    end

    if all_in
        push!(uref_coll, refs[2])
    end
end

# "Set" (vs "add") assumes U_sub is empty initially... May fail if not respected
# function set_inputs_by_xref_ysub!(U_sub::SubSetHash, sym_model::SymbolicModelHash, x_ref, Y_sub)
#     add_inputs_by_xref_ysub!(U_sub, sym_model, x_ref, Y_sub)
#     U_sub.issorted = true
#     U_sub.isunique = true
# end

# function add_inputs_images_by_xref!(U_sub, Y_sub, sym_model::SymbolicModelHash, x_ref)
#     ensure_sorted!(sym_model)
#     ensure_unique!(sym_model)
#     idx_iter = searchsorted(sym_model.elems, (x_ref,), by = x -> x[1])
#     uref_prev = U_sub.grid_space.overflow_ref
#
#     for idx in idx_iter
#         refs = sym_model.elems[idx]
#         if refs[2] != uref_prev
#             add_to_subset_by_ref!(U_sub, refs[2])
#             uref_prev = refs[2]
#         end
#         add_to_subset_by_ref!(Y_sub, refs[3])
#     end
# end

function add_inputs_images_by_xref!(uref_coll, yref_coll, sym_model::SymbolicModelHash, x_ref)
    ensure_sorted!(sym_model)
    ensure_unique!(sym_model)
    idx_iter = searchsorted(sym_model.elems, (x_ref,), by = x -> x[1])
    uref_prev = sym_model.U_grid.overflow_ref

    for idx in idx_iter
        refs = sym_model.elems[idx]
        if refs[2] != uref_prev
            push!(uref_coll, refs[2])
            uref_prev = refs[2]
        end
        push!(yref_coll, refs[3])
    end
end

function get_symmodel_size(sym_model::SymbolicModelHash)
    ensure_unique!(sym_model)
    return length(sym_model.elems)
end

function sizehint_symmodel!(sym_model::SymbolicModelHash, size_max)
    sizehint!(sym_model.elems, size_max)
end
