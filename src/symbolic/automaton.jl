abstract type AbstractAutomatonList{N, M} <: HybridSystems.AbstractAutomaton end

# === Required Interface ===
function get_n_state(autom::AbstractAutomatonList{N, M}) where {N, M} end
function get_n_input(autom::AbstractAutomatonList{N, M}) where {N, M} end
function enum_transitions(autom::AbstractAutomatonList{N, M}) where {N, M} end
function add_transition!(
    autom::AbstractAutomatonList{N, M},
    source::Int,
    target::Int,
    symbol::Int,
) where {N, M} end
function pre(autom::AbstractAutomatonList{N, M}, target::Int) where {N, M} end
function post(autom::AbstractAutomatonList{N, M}, source::Int, symbol::Int) where {N, M} end
function empty!(autom::AbstractAutomatonList{N, M}) where {N, M} end
function add_state!(autom::AbstractAutomatonList{N, M}) where {N, M} end

# === Common Default Implementations ===
enum_states(autom::AbstractAutomatonList) = 1:get_n_state(autom)
enum_inputs(autom::AbstractAutomatonList) = 1:get_n_input(autom)

function HybridSystems.ntransitions(autom::AbstractAutomatonList{N, M}) where {N, M}
    return length(enum_transitions(autom))
end

function add_transitions!(autom::AbstractAutomatonList{N, M}, translist) where {N, M}
    for (q′, q, u) in translist
        add_transition!(autom, q, q′, u)
    end
end

function is_deterministic(autom::AbstractAutomatonList{N, M}) where {N, M}
    seen = Dict{Tuple{Int, Int}, Int}()
    for (q′, q, u) in enum_transitions(autom)
        key = (q, u)
        seen[key] = get(seen, key, 0) + 1
        if seen[key] > 1
            return false
        end
    end
    return true
end

# === Implementation: SortedAutomatonList ===

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

function post(a::SortedAutomatonList, source::Int, symbol::Int)
    targets = Int[]
    UT.fix_and_eliminate_tail!(targets, a.transitions, (source, symbol))
    return targets
end

Base.empty!(a::SortedAutomatonList) = empty!(a.transitions)
function HybridSystems.add_state!(a::SortedAutomatonList)
    a.nstates += 1
    return a.nstates
end

# === Implementation: IndexedAutomatonList ===

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
