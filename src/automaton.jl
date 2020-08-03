abstract type Automaton end

mutable struct AutomatonList <: Automaton
    nstates::Int
    nsymbols::Int
    transitions::Vector{Tuple{Int, Int, Int}}
    issorted::Bool
end

function NewAutomatonList(nstates, nsymbols)
    transitions = Tuple{Int, Int, Int}[]
    return AutomatonList(nstates, nsymbols, transitions, true)
end

function ensure_sorted!(autom::AutomatonList)
    if !autom.issorted
        # display("autom not sorted")
        sort!(autom.transitions)
        autom.issorted = true
    end
end

function get_ntrans(autom::AutomatonList)
    return length(autom.transitions)
end

# Do not check that source, symbol, target are "inbounds"
# Assumes not add twice same transition...
function add_transition!(autom::AutomatonList, source, symbol, target)
    push!(autom.transitions, (target, source, symbol))
    autom.issorted = false
end

function Base.empty!(autom::AutomatonList)
    empty!(autom.transitions)
    autom.issorted = true
end

function compute_post!(targetlist, autom::AutomatonList, source, symbol)
    ensure_sorted!(autom)
    check = x -> x[2:3] == (source, symbol)
    for trans in Iterators.filter(check, autom.transitions)
        push!(targetlist, trans[1])
    end
end

function compute_pre!(soursymblist, autom::AutomatonList, target)
    ensure_sorted!(autom)
    idxlist = searchsorted(autom.transitions, (target,), by = x -> x[1])
    for idx in idxlist
        push!(soursymblist, autom.transitions[idx][2:3])
    end
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
