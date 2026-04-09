mutable struct IndexedAutomatonList <: AbstractAutomatonList{3, 3}
    nstates::Int
    nsymbols::Int
    transitions::Vector{NTuple{3, Int}}
    postmap::Dict{Tuple{Int, Int}, Vector{Int}}
    premap::Dict{Int, Vector{Tuple{Int, Int}}}
end

function NewIndexedAutomatonList(nstates, nsymbols)
    return IndexedAutomatonList(nstates, nsymbols, NTuple{3, Int}[], Dict(), Dict())
end

get_n_state(a::IndexedAutomatonList) = a.nstates
get_n_input(a::IndexedAutomatonList) = a.nsymbols
enum_transitions(a::IndexedAutomatonList) = a.transitions

function add_transition!(a::IndexedAutomatonList, q::Int, q′::Int, u::Int)
    push!(a.transitions, (q′, q, u))
    push!(get!(a.postmap, (q, u), Int[]), q′)
    return push!(get!(a.premap, q′, Tuple{Int, Int}[]), (q, u))
end

function add_transitions!(autom::IndexedAutomatonList, translist)
    append!(autom.transitions, translist)
    for (q′, q, u) in translist
        push!(get!(autom.postmap, (q, u), Int[]), q′)
        push!(get!(autom.premap, q′, Tuple{Int, Int}[]), (q, u))
    end
end

post(a::IndexedAutomatonList, q::Int, u::Int) = get(a.postmap, (q, u), Int[])
pre(a::IndexedAutomatonList, q′::Int) = get(a.premap, q′, Tuple{Int, Int}[])

Base.empty!(a::IndexedAutomatonList) = begin
    empty!(a.transitions)
    empty!(a.postmap)
    empty!(a.premap)
end
function HybridSystems.add_state!(a::IndexedAutomatonList)
    a.nstates += 1
    return a.nstates
end
