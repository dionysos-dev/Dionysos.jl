"""
    FastIndexedAutomatonList <: AbstractAutomatonList

Automaton with dense vector indices: `postmap` is keyed by the flattened
`(state, symbol)` pair and `premap` by target state, so `post`/`pre` are direct
array lookups. Fastest to query and densest in memory; best for large fixed
automata. Call [`finalize!`](@ref) after bulk insertion to deduplicate.
"""
mutable struct FastIndexedAutomatonList <: AbstractAutomatonList
    nstates::Int
    nsymbols::Int
    transitions::Vector{NTuple{3, Int}}
    postmap::Vector{Vector{Int}}
    premap::Vector{Vector{Tuple{Int, Int}}}
end

function FastIndexedAutomatonList(nstates::Int, nsymbols::Int)
    return FastIndexedAutomatonList(
        nstates,
        nsymbols,
        NTuple{3, Int}[],
        [Int[] for _ in 1:(nstates * nsymbols)],
        [Tuple{Int, Int}[] for _ in 1:nstates],
    )
end

_pair_id(a::FastIndexedAutomatonList, q::Int, u::Int) = (q - 1) * a.nsymbols + u

get_n_state(a::FastIndexedAutomatonList) = a.nstates
get_n_input(a::FastIndexedAutomatonList) = a.nsymbols
enum_transitions(a::FastIndexedAutomatonList) = a.transitions

function add_transition!(a::FastIndexedAutomatonList, q::Int, q′::Int, u::Int)
    push!(a.transitions, (q′, q, u))
    push!(a.postmap[_pair_id(a, q, u)], q′)
    push!(a.premap[q′], (q, u))
    return a
end

function add_transitions!(a::FastIndexedAutomatonList, translist)
    sizehint!(a.transitions, length(a.transitions) + length(translist))
    append!(a.transitions, translist)

    for (q′, q, u) in translist
        push!(a.postmap[_pair_id(a, q, u)], q′)
        push!(a.premap[q′], (q, u))
    end

    return a
end

post(a::FastIndexedAutomatonList, q::Int, u::Int) = a.postmap[_pair_id(a, q, u)]

pre(a::FastIndexedAutomatonList, q′::Int) = a.premap[q′]

function Base.empty!(a::FastIndexedAutomatonList)
    empty!(a.transitions)
    for v in a.postmap
        empty!(v)
    end
    for v in a.premap
        empty!(v)
    end
    return a
end

function HybridSystems.add_state!(a::FastIndexedAutomatonList)
    a.nstates += 1

    for _ in 1:a.nsymbols
        push!(a.postmap, Int[])
    end

    push!(a.premap, Tuple{Int, Int}[])

    return a.nstates
end

function finalize!(autom::FastIndexedAutomatonList)
    for v in autom.postmap
        length(v) > 1 && unique!(v)
    end

    for v in autom.premap
        length(v) > 1 && unique!(v)
    end

    length(autom.transitions) > 1 && unique!(autom.transitions)

    return autom
end
