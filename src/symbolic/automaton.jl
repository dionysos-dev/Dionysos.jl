mutable struct AutomatonList <: HybridSystems.AbstractAutomaton
    nstates::Int
    nsymbols::Int
    transitions::UT.SortedTupleSet{3,Int}
end

function NewAutomatonList(nstates, nsymbols)
    transitions = UT.SortedTupleSet{3,Int}()
    return AutomatonList(nstates, nsymbols, transitions)
end

function HybridSystems.ntransitions(autom::AutomatonList)
    return length(autom.transitions)
end

# In add_trans and add_translist:
# Do not check that source, symbol, target are "inbounds"
# Assumes not add twice same transition...
function HybridSystems.add_transition!(autom::AutomatonList, source, target, symbol)
    UT.push_new!(autom.transitions, (target, source, symbol))
end

# translist is an iterable of Tuple{Int,Int,Int}
function add_transitions!(autom::AutomatonList, translist)
    UT.append_new!(autom.transitions, translist)
end
Base.empty!(autom::AutomatonList) = empty!(autom.transitions)

function compute_post!(targetlist, autom::AutomatonList, source, symbol)
    UT.fix_and_eliminate_tail!(targetlist, autom.transitions, (source, symbol))
end

function pre(autom::AutomatonList, target)
    return UT.fix_and_eliminate_first(autom.transitions, target)
end

# function add_inputs_images_by_xref!(uref_coll, yref_coll, autom::AutomatonList, x_ref)
#     ensure_sorted!(autom)
#     ensure_unique!(autom)
#     idx_iter = searchsorted(autom.transitions, (x_ref,), by = x -> x[1])
#     uref_prev = autom.Ugrid.overflow_ref
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
