# ============================================================
# Lifted quotient data structures
# ============================================================

mutable struct PCAbstractState{S,U}
    id::Int
    node::U
    set::S
    obs::Int
    slice::Int
    next::Vector{Tuple{Int,Int}}   # (mode, target_state_id)
end

mutable struct PCBisimulationQuotient{S,U}
    states::Dict{Int,PCAbstractState{S,U}}
    part_ids::Dict{U,Vector{Int}}   # state ids per node partition
    next_id::Int
end

function PCBisimulationQuotient{S,U}() where {S,U}
    return PCBisimulationQuotient{S,U}(
        Dict{Int,PCAbstractState{S,U}}(),
        Dict{U,Vector{Int}}(),
        1,
    )
end

function add_state!(
    T::PCBisimulationQuotient{S,U},
    node::U,
    set::S,
    obs::Int,
    slice::Int,
) where {S,U}
    qid = T.next_id
    T.next_id += 1
    q = PCAbstractState{S,U}(qid, node, set, obs, slice, Tuple{Int,Int}[])
    T.states[qid] = q
    get!(T.part_ids, node, Int[])
    push!(T.part_ids[node], qid)
    return qid
end

function remove_state!(T::PCBisimulationQuotient, qid::Int)
    q = T.states[qid]
    if haskey(T.part_ids, q.node)
        filter!(x -> x != qid, T.part_ids[q.node])
    end
    delete!(T.states, qid)
    return nothing
end

function add_transition!(T::PCBisimulationQuotient, src::Int, mode::Int, dst::Int)
    tr = (mode, dst)
    if tr ∉ T.states[src].next
        push!(T.states[src].next, tr)
    end
    return nothing
end