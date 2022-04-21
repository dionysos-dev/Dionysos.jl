mutable struct AutomatonList{S<:AbstractSet{NTuple{3,Int}}} <: HybridSystems.AbstractAutomaton
    nstates::Int
    nsymbols::Int
    transitions::S
end

function AutomatonList{S}(nstates, nsymbols) where {S}
    transitions = S()
    return AutomatonList(nstates, nsymbols, transitions)
end
NewAutomatonList(nstates, nsymbols) = AutomatonList{SortedTupleSet{3,Int}}(nstates, nsymbols)

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
append_new!(s::Set, translist) = union!(s, translist)
function add_transitions!(autom::AutomatonList, translist)
    UT.append_new!(autom.transitions, translist)
end

Base.empty!(autom::AutomatonList) = empty!(autom.transitions)

function compute_post!(targetlist, autom::AutomatonList, source, symbol)
    UT.fix_and_eliminate_tail!(targetlist, autom.transitions, (source, symbol))
end

#function compute_available!(targetlist, autom::AutomatonList, source)
   # fix_and_eliminate_tail!(targetlist, autom.transitions, (source, symbol))
#end

function pre(autom::AutomatonList, target)
    return UT.fix_and_eliminate_first(autom.transitions, target)
end

function HybridSystems.add_state!(autom::AutomatonList)
    autom.nstates += 1
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
