function NewAutomatonList(nstates, nsymbols)
    return AutomatonList(nstates, nsymbols, Tuple{CellRef, CellRef, CellRef}[], true, true)
end

function ensure_sorted!(autom::AutomatonList)
    if !autom.issorted
        # display("autom not sorted")
        sort!(autom.transitions)
        autom.issorted = true
    end
end

# Do not check that x, u, y are "inbounds"
# Assume not add twice same transition...
function add_transition!(autom::AutomatonList, source, symbol, target)
    push!(autom.transitions, (target, source, symbol))
    autom.issorted = false
end

function empty!(autom::AutomatonList)
    empty!(autom.transitions)
    autom.issorted = true
end

function compute_post!(targetlist, autom::AutomatonList, source, symbol)
    ensure_sorted!(autom)
    for trans in Iterators.filter(x -> x[2:3] == (source, symbol), autom.transitions)
        push!(targetlist, trans[1])
    end
end

function compute_pre!(sourcesymbollist, autom::AutomatonList, target)
    ensure_sorted!(autom)
    idxlist = searchsorted(autom.transitions, (target,), by = x -> x[1])
    for idx in idxlist
        push!(sourcesymbollist, autom.transitions[idx][2:3])
    end
end

# function add_inputs_images_by_xref!(uref_coll, yref_coll, autom::AutomatonList, x_ref)
#     ensure_sorted!(autom)
#     ensure_unique!(autom)
#     idx_iter = searchsorted(autom.transitions, (x_ref,), by = x -> x[1])
#     uref_prev = autom.U_grid.overflow_ref
#
#     for idx in idx_iter
#         refs = autom.transitions[idx]
#         if refs[2] != uref_prev
#             push!(uref_coll, refs[2])
#             uref_prev = refs[2]
#         end
#         push!(yref_coll, refs[3])
#     end
# end
#
# function get_symmodel_size(autom::AutomatonList)
#     ensure_unique!(autom)
#     return length(autom.transitions)
# end
#
# function sizehint_symmodel!(autom::AutomatonList, size_max)
#     sizehint!(autom.transitions, size_max)
# end
