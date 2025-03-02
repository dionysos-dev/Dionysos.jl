mutable struct AutomatonList{S <: AbstractSet{NTuple{3, Int}}} <:
               HybridSystems.AbstractAutomaton
    nstates::Int
    nsymbols::Int
    transitions::S
end

function AutomatonList{S}(nstates, nsymbols) where {S}
    transitions = S()
    return AutomatonList(nstates, nsymbols, transitions)
end

NewAutomatonList(nstates, nsymbols) =
    AutomatonList{UT.SortedTupleSet{3, NTuple{3, Int}}}(nstates, nsymbols)

function HybridSystems.ntransitions(autom::AutomatonList)
    return length(autom.transitions)
end

# Assumes not add twice same transition...
function HybridSystems.add_transition!(autom::AutomatonList, source, target, symbol)
    return UT.push_new!(autom.transitions, (target, source, symbol))
end

function add_transitions!(autom::AutomatonList, translist)
    return UT.append_new!(autom.transitions, translist)
end

Base.empty!(autom::AutomatonList) = empty!(autom.transitions)

function compute_post!(targetlist, autom::AutomatonList, source, symbol)
    return UT.fix_and_eliminate_tail!(targetlist, autom.transitions, (source, symbol))
end

function pre(autom::AutomatonList, target)
    return UT.fix_and_eliminate_first(autom.transitions, target)
end

function HybridSystems.add_state!(autom::AutomatonList)
    return autom.nstates += 1
end

function is_deterministic(autom::AutomatonList)
    seen = Dict{Tuple{Int, Int}, Int}()  # Stores (source, symbol) → number of targets

    for (target, source, symbol) in UT.get_data(autom.transitions)
        key = (source, symbol)
        seen[key] = get(seen, key, 0) + 1
        if seen[key] > 1  # If we find more than one target for (source, symbol), it's non-deterministic
            return false
        end
    end
    return true
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
