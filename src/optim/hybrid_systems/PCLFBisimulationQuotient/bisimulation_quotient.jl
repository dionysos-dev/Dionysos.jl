# ============================================================
# Lifted quotient data structures
# ============================================================

mutable struct PCAbstractState{S, U}
    id::Int
    node::U
    set::S
    obs::Int
    slice::Int
    next::Vector{Tuple{Int, Int}}   # (mode, target_state_id)
end

mutable struct PCBisimulationQuotient{S, U}
    states::Dict{Int, PCAbstractState{S, U}}
    part_ids::Dict{U, Vector{Int}}
    next_id::Int
    slices::Dict{U, Vector{S}}
end

function PCBisimulationQuotient{S, U}(slices::Dict{U, Vector{S}}) where {S, U}
    return PCBisimulationQuotient{S, U}(
        Dict{Int, PCAbstractState{S, U}}(),
        Dict{U, Vector{Int}}(),
        1,
        slices,
    )
end

function add_state!(
    T::PCBisimulationQuotient{S, U},
    node::U,
    set::S,
    obs::Int,
    slice::Int,
) where {S, U}
    qid = T.next_id
    T.next_id += 1
    q = PCAbstractState{S, U}(qid, node, set, obs, slice, Tuple{Int, Int}[])
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

function semilinear_by_node(T::PCBisimulationQuotient{S, U}, state_ids) where {S, U}
    parts_by_node = Dict{U, Vector{Poly}}()

    for qid in state_ids
        haskey(T.states, qid) || continue
        q = T.states[qid]
        get!(parts_by_node, q.node, Poly[])
        append!(parts_by_node[q.node], q.set.array)
    end

    return Dict(
        nd => UT.normalize_semilinear(UT.semilinear_set(parts)) for
        (nd, parts) in parts_by_node
    )
end

# Bounding boxes of a batch of parts, laid out column-wise: one pair of matrices per batch
# rather than a `(low, high)` vector pair per part, so each box is contiguous and the whole
# batch costs a single allocation. Dimension is read off the parts, so this stays
# independent of the state-space dimension.
function _part_boxes(ps)
    n = isempty(ps) ? 0 : LazySets.dim(first(ps))
    lo = Matrix{Float64}(undef, n, length(ps))
    hi = Matrix{Float64}(undef, n, length(ps))
    for (j, P) in enumerate(ps)
        H = LazySets.box_approximation(P)
        l, h = LazySets.low(H), LazySets.high(H)
        @inbounds for i in 1:n
            lo[i, j] = l[i]
            hi[i, j] = h[i]
        end
    end
    return lo, hi
end

@inline function _boxes_overlap(lo_a, hi_a, ja, lo_b, hi_b, jb)
    @inbounds for i in axes(lo_a, 1)
        (lo_a[i, ja] <= hi_b[i, jb] && lo_b[i, jb] <= hi_a[i, ja]) || return false
    end
    return true
end

# Nonempty parts of the pairwise intersection, box-screened. Kept a separate function on
# purpose: the caller reassigns the running intersection and its boxes as it chains through
# the nodes, which leaves the compiler boxing those locals, and the box lookup would then
# land on every one of the |parts|·|other| screens. Passing them as arguments restores
# concrete types for the sweep.
function _intersect_all(parts, lo_a, hi_a, other, lo_b, hi_b)
    out = similar(parts, 0)
    for ja in eachindex(parts), jb in eachindex(other)
        _boxes_overlap(lo_a, hi_a, ja, lo_b, hi_b, jb) || continue
        I = UT.poly_intersection(parts[ja], other[jb])
        if I isa LazySets.HPolytope && !isempty(I)
            push!(out, I)
        end
    end
    return out
end

"""
    get_volume(T::PCBisimulationQuotient, state_ids; backend, atol = 1e-6)

Volume of the region covered by `state_ids`, counting the overlap between nodes once.

Parts sharing a node are pairwise disjoint — they are equivalence classes of the quotient —
so the union volume is inclusion-exclusion over the nodes rather than a running set
difference. That matters twice over. `set_difference_decompose` emits one piece per
constraint of its subtrahend and feeds those pieces into the next subtraction, so the piece
count compounds; an intersection yields at most one piece and never feeds back. And a
complement halfspace is inset by `atol`, eroding a sliver at every cut, so thousands of cuts
accumulate a one-sided deficit that grows with the quotient. Intersections form no
complement, so the result no longer depends on `atol` at all.

`atol` is kept for interface compatibility and does not affect the value.
"""
function get_volume(
    T::PCBisimulationQuotient,
    state_ids;
    backend = nothing,
    atol::Float64 = 1e-6,
)
    backend === nothing &&
        error("No polyhedral backend provided. Example: backend = CDDLib.Library().")

    S_by_node = semilinear_by_node(T, state_ids)
    isempty(S_by_node) && return 0.0

    node_parts = [Snode.array for Snode in values(S_by_node)]
    sum_volume(ps) =
        isempty(ps) ? 0.0 : sum(LazySets.volume(P; backend = backend) for P in ps)

    total = sum(sum_volume(ps) for ps in node_parts)
    M = length(node_parts)
    M == 1 && return total

    node_boxes = map(_part_boxes, node_parts)

    # Subtract pairwise overlaps, add back triples, and so on. The higher-order terms die
    # out quickly because the running intersection empties, so the 2^M bound on the number
    # of subsets is not approached in practice.
    for mask in 1:((1 << M) - 1)
        nodes = [i for i in 1:M if (mask >> (i - 1)) & 1 == 1]
        length(nodes) < 2 && continue

        parts = node_parts[nodes[1]]
        lo_a, hi_a = node_boxes[nodes[1]]
        for (step, k) in enumerate(@view nodes[2:end])
            other = node_parts[k]
            lo_b, hi_b = node_boxes[k]
            parts = _intersect_all(parts, lo_a, hi_a, other, lo_b, hi_b)
            isempty(parts) && break
            # Boxes of the running intersection are only worth building when a further node
            # still has to be screened against them.
            if step < length(nodes) - 1
                lo_a, hi_a = _part_boxes(parts)
            end
        end

        total += (isodd(length(nodes)) ? 1.0 : -1.0) * sum_volume(parts)
    end

    return total
end
