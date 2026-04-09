mutable struct SortedAutomatonList{S <: AbstractSet{NTuple{3, Int}}} <:
               AbstractAutomatonList{3, 3}
    nstates::Int
    nsymbols::Int
    transitions::S
end

function SortedAutomatonList{S}(nstates, nsymbols) where {S}
    return SortedAutomatonList(nstates, nsymbols, S())
end

function NewSortedAutomatonList(nstates, nsymbols)
    return SortedAutomatonList{UT.SortedTupleSet{3, NTuple{3, Int}}}(nstates, nsymbols)
end

get_n_state(a::SortedAutomatonList) = a.nstates
get_n_input(a::SortedAutomatonList) = a.nsymbols
enum_transitions(a::SortedAutomatonList) = UT.get_data(a.transitions)
add_transition!(a::SortedAutomatonList, q::Int, q′::Int, u::Int) =
    UT.push_new!(a.transitions, (q′, q, u))
add_transitions!(autom::SortedAutomatonList, translist) =
    UT.append_new!(autom.transitions, translist)
pre(a::SortedAutomatonList, target::Int) = UT.fix_and_eliminate_first(a.transitions, target)

compute_post!(targetlist, a::SortedAutomatonList, source, symbol) =
    UT.fix_and_eliminate_tail!(targetlist, a.transitions, (source, symbol))

function post(a::SortedAutomatonList, source::Int, symbol::Int)
    targets = Int[]
    UT.fix_and_eliminate_tail!(targets, a.transitions, (source, symbol))
    return targets
end

Base.empty!(a::SortedAutomatonList) = Base.empty!(a.transitions)
function HybridSystems.add_state!(a::SortedAutomatonList)
    a.nstates += 1
    return a.nstates
end
