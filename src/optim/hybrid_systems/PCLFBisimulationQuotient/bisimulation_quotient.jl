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

function semilinear_by_node(T::PCBisimulationQuotient, state_ids)
    parts_by_node = Dict{Any, Vector{Poly}}()

    for qid in state_ids
        haskey(T.states, qid) || continue
        q = T.states[qid]
        get!(parts_by_node, q.node, Poly[])
        append!(parts_by_node[q.node], q.set.parts)
    end

    return Dict(
        nd => UT.normalize_semilinear(UT.SemiLinearSet(parts)) for
        (nd, parts) in parts_by_node
    )
end

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

    if length(S_by_node) == 1
        return sum(
            UT.get_volume(P; backend = backend) for P in first(values(S_by_node)).parts
        )
    end

    total = 0.0
    accumulated_parts = Poly[]

    for Snode in values(S_by_node)
        current = copy(Snode.parts)

        for Q in accumulated_parts
            isempty(current) && break

            new_current = Poly[]
            for P in current
                I = UT.set_intersection(P, Q)
                if UT.is_nonempty_set(I)
                    append!(new_current, UT.set_difference_decompose(P, Q; atol = atol))
                else
                    push!(new_current, P)
                end
            end
            current = new_current
        end

        total += sum(UT.get_volume(P; backend = backend) for P in current)
        append!(accumulated_parts, Snode.parts)
    end

    return total
end

# ------------------------------------------------------------
# Single abstract state
# ------------------------------------------------------------
@recipe function f(q::PCAbstractState; color = :blue, fillalpha = 0.25, show_label = true)
    linecolor := color
    fillcolor := color
    fillalpha := fillalpha
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

            for (j, P) in enumerate(q.set.parts)
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

            for (j, P) in enumerate(S.parts)
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
    return length(S.parts)
end

function num_faces(P::UT.Poly)
    Q = UT.clean_poly(copy(P))
    return length(LazySets.constraints_list(Q))
end

function num_faces(S::UT.SemiLinearSet)
    total = 0
    for P in S.parts
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
