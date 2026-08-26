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

# ------------------------------------------------------------
# Single abstract state
# ------------------------------------------------------------
# The styling is declared with `-->` rather than taken as keyword arguments: Plots rewrites its
# own attribute aliases (`color` becomes `seriescolor`) before a recipe runs, so a keyword named
# `color` here could never be read. As defaults they behave the ordinary way — the caller wins.
@recipe function f(q::PCAbstractState; show_label = true)
    seriescolor --> :blue
    fillalpha --> 0.25
    label := show_label ? "q$(q.id) [node=$(q.node), slice=$(q.slice), obs=$(q.obs)]" : ""
    return q.set
end

# ------------------------------------------------------------
# Full quotient recipe
# ------------------------------------------------------------
@recipe function f(
    T::PCBisimulationQuotient;
    what = :states,
    node = nothing,
    slice = nothing,
    obs = nothing,
    mode = nothing,
    state_ids = nothing,
    by = :state,
    fillalpha = 0.25,
    linewidth = 1.5,
    seriesalpha = 0.9,
    show_labels = false,
    show_contours = false,
    user_color = nothing,
)
    palette = [
        :red,
        :blue,
        :green,
        :orange,
        :purple,
        :brown,
        :pink,
        :cyan,
        :magenta,
        :olive,
        :gold,
        :coral,
        :turquoise,
        :navy,
        :darkgreen,
        :darkred,
    ]

    local_linealpha = show_contours ? 1.0 : 0.0
    state_id_set = isnothing(state_ids) ? nothing : Set(state_ids)

    if what == :states
        qlist = [
            q for q in values(T.states) if (isnothing(node) || q.node == node) &&
                (isnothing(slice) || q.slice == slice) &&
                (isnothing(obs) || q.obs == obs) &&
                (isnothing(state_id_set) || q.id in state_id_set)
        ]

        sort!(qlist; by = q -> q.id)

        for (k, q) in enumerate(qlist)
            c = if !isnothing(user_color)
                user_color
            elseif by == :slice
                palette[mod1(q.slice, length(palette))]
            elseif by == :obs
                palette[mod1(q.obs + 2, length(palette))]
            elseif by == :node
                palette[mod1(abs(hash(q.node)), length(palette))]
            else
                palette[mod1(k, length(palette))]
            end

            for (j, P) in enumerate(q.set.array)
                @series begin
                    seriestype := :shape
                    fillcolor := c
                    linecolor := c
                    fillalpha := fillalpha
                    linealpha := local_linealpha
                    linewidth := linewidth
                    seriesalpha := seriesalpha
                    label := (show_labels && j == 1) ? "q$(q.id)" : ""
                    P
                end
            end
        end
    end

    if what == :slices
        seen = Set{Tuple{Any, Int}}()

        groups = Tuple{Any, Int, UT.SemiLinearSet}[]
        for (nd, slice_list) in T.slices
            if !isnothing(node) && nd != node
                continue
            end
            for (i, S) in enumerate(slice_list)
                if !isnothing(slice) && i != slice
                    continue
                end
                push!(groups, (nd, i, S))
            end
        end

        sort!(groups; by = x -> x[2])

        for (nd, i, S) in groups
            c = palette[mod1(i, length(palette))]
            key = (nd, i)

            for (j, P) in enumerate(S.array)
                @series begin
                    seriestype := :shape
                    fillcolor := c
                    linecolor := c
                    fillalpha := fillalpha
                    linealpha := local_linealpha
                    linewidth := linewidth
                    seriesalpha := seriesalpha
                    label :=
                        (show_labels && !(key in seen) && j == 1) ? "node=$nd, slice=$i" :
                        ""
                    P
                end
            end

            push!(seen, key)
        end
    end
end

# ============================================================
# Basic statistics on the bisimulation quotient
# ============================================================

num_nodes(T::PCBisimulationQuotient) = length(T.slices)
function num_slices(T::PCBisimulationQuotient)
    isempty(T.slices) && return 0
    return length(first(values(T.slices)))
end

num_states(T::PCBisimulationQuotient) = length(T.states)
num_transitions(T::PCBisimulationQuotient) = sum(length(q.next) for q in values(T.states))

function states_by_obs(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for q in values(T.states)
        d[q.obs] = get(d, q.obs, 0) + 1
    end
    return sort(collect(d); by = first) |> Dict
end

function states_by_slice(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for q in values(T.states)
        d[q.slice] = get(d, q.slice, 0) + 1
    end
    return sort(collect(d); by = first) |> Dict
end

function states_by_node(T::PCBisimulationQuotient)
    d = Dict{typeof(first(keys(T.part_ids))), Int}()
    for q in values(T.states)
        d[q.node] = get(d, q.node, 0) + 1
    end
    return d
end

function transitions_by_mode(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for q in values(T.states)
        for (m, _) in q.next
            d[m] = get(d, m, 0) + 1
        end
    end
    return sort(collect(d); by = first) |> Dict
end

function outgoing_degree_by_state(T::PCBisimulationQuotient)
    d = Dict{Int, Int}()
    for (qid, q) in T.states
        d[qid] = length(q.next)
    end
    return d
end

function outgoing_degree_stats(T::PCBisimulationQuotient)
    degs = [length(q.next) for q in values(T.states)]
    isempty(degs) && return Dict(:min => 0, :max => 0, :mean => 0.0, :median => 0.0)
    return Dict(
        :min => minimum(degs),
        :max => maximum(degs),
        :mean => sum(degs) / length(degs),
        :median => sort(degs)[cld(length(degs), 2)],
    )
end

function states_in_slice(T::PCBisimulationQuotient, slice::Int; states = values(T.states))
    return [q for q in states if q.slice == slice]
end

function states_in_obs(T::PCBisimulationQuotient, obs::Int; states = values(T.states))
    return [q for q in states if q.obs == obs]
end

function states_in_node(T::PCBisimulationQuotient, node; states = values(T.states))
    return [q for q in states if q.node == node]
end

function state_ids_in_slice(
    T::PCBisimulationQuotient,
    slice::Int;
    state_ids = keys(T.states),
)
    return [
        qid for qid in state_ids if haskey(T.states, qid) && T.states[qid].slice == slice
    ]
end

function state_ids_in_obs(T::PCBisimulationQuotient, obs::Int; state_ids = keys(T.states))
    return [qid for qid in state_ids if haskey(T.states, qid) && T.states[qid].obs == obs]
end

function state_ids_in_node(T::PCBisimulationQuotient, node; state_ids = keys(T.states))
    return [qid for qid in state_ids if haskey(T.states, qid) && T.states[qid].node == node]
end

function transitions_from_mode(T::PCBisimulationQuotient, mode::Int)
    out = Vector{Tuple{Int, Int}}()
    for (qid, q) in T.states
        for (m, dst) in q.next
            if m == mode
                push!(out, (qid, dst))
            end
        end
    end
    return out
end

function reachable_target_states(T::PCBisimulationQuotient)
    targets = Set{Int}()
    for q in values(T.states)
        for (_, dst) in q.next
            push!(targets, dst)
        end
    end
    return targets
end

function deadend_states(T::PCBisimulationQuotient)
    return [q.id for q in values(T.states) if isempty(q.next)]
end

function self_loop_count(T::PCBisimulationQuotient)
    c = 0
    for (qid, q) in T.states
        for (_, dst) in q.next
            c += (dst == qid)
        end
    end
    return c
end

function num_parts(S::UT.SemiLinearSet)
    return length(S.array)
end

function num_faces(P::Poly)
    Q = UT.clean_poly(copy(P))
    return length(LazySets.constraints_list(Q))
end

function num_faces(S::UT.SemiLinearSet)
    total = 0
    for P in S.array
        try
            total += num_faces(P)
        catch
            continue
        end
    end
    return total
end

function cell_complexities(T::PCBisimulationQuotient)
    n_parts = Int[]
    n_faces = Int[]

    for q in values(T.states)
        S = q.set
        push!(n_parts, num_parts(S))
        push!(n_faces, num_faces(S))
    end

    return n_parts, n_faces
end

function bisimulation_stats(T::PCBisimulationQuotient)
    return Dict(
        :num_states => num_states(T),
        :num_transitions => num_transitions(T),
        :states_by_obs => states_by_obs(T),
        :states_by_slice => states_by_slice(T),
        :states_by_node => states_by_node(T),
        :transitions_by_mode => transitions_by_mode(T),
        :outgoing_degree_stats => outgoing_degree_stats(T),
        :deadend_states => deadend_states(T),
        :num_deadend_states => length(deadend_states(T)),
        :self_loop_count => self_loop_count(T),
    )
end

function print_bisimulation_stats(T::PCBisimulationQuotient)
    println("Bisimulation quotient statistics")
    println("--------------------------------")
    println("Number of nodes        : ", num_nodes(T))
    println("Number of slices       : ", num_slices(T))
    println("Number of states       : ", num_states(T))
    println("Number of transitions  : ", num_transitions(T))

    println("States by observation  : ", sort(collect(states_by_obs(T)); by = first))
    println("States by slice        : ", sort(collect(states_by_slice(T)); by = first))
    println("States by node         : ", sort(collect(states_by_node(T)); by = first))
    println("Transitions by mode    : ", sort(collect(transitions_by_mode(T)); by = first))

    println("Outgoing degree stats  : ", outgoing_degree_stats(T))
    println("Deadend states         : ", length(deadend_states(T)))
    println("Self-loops             : ", self_loop_count(T))
    return nothing
end
