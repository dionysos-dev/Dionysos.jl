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
    part_ids::Dict{U, Vector{Int}}  # state ids per node partition
    next_id::Int
    slices::Dict{U, Vector{Vector{S}}}
    obs_partition::Vector{Tuple{S, Int}}
end

function PCBisimulationQuotient{S, U}(
    slices::Dict{U, Vector{Vector{S}}},
    obs_partition::Vector{Tuple{S, Int}},
) where {S, U}
    return PCBisimulationQuotient{S, U}(
        Dict{Int, PCAbstractState{S, U}}(),
        Dict{U, Vector{Int}}(),
        1,
        slices,
        obs_partition,
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
    by = :slice,
    fillalpha = 0.25,
    linewidth = 1.5,
    seriesalpha = 0.9,
    show_labels = false,
)
    palette =
        [:red, :blue, :green, :orange, :purple, :brown, :pink, :cyan, :magenta, :olive]

    # --------------------------------------------------------
    # Plot abstract states
    # --------------------------------------------------------
    if what == :states
        for q in values(T.states)
            if !isnothing(node) && q.node != node
                continue
            end
            if !isnothing(slice) && q.slice != slice
                continue
            end
            if !isnothing(obs) && q.obs != obs
                continue
            end

            c = if by == :slice
                palette[mod1(q.slice, length(palette))]
            elseif by == :obs
                palette[mod1(q.obs + 2, length(palette))]
            elseif by == :node
                palette[mod1(abs(hash(q.node)), length(palette))]
            else
                :blue
            end

            @series begin
                linecolor := c
                fillcolor := c
                fillalpha := fillalpha
                label := show_labels ? "q$(q.id)" : ""
                q.set
            end
        end
    end

    # --------------------------------------------------------
    # Plot slices stored in the quotient
    # --------------------------------------------------------
    if what == :slices
        seen = Set{Tuple{Any, Int}}()

        for (nd, slice_list) in T.slices
            if !isnothing(node) && nd != node
                continue
            end

            for (i, polys) in enumerate(slice_list)
                if !isnothing(slice) && i != slice
                    continue
                end

                c = palette[mod1(i, length(palette))]
                key = (nd, i)

                for P in polys
                    @series begin
                        linecolor := c
                        fillcolor := c
                        fillalpha := fillalpha
                        label := (show_labels && !(key in seen)) ? "node=$nd, slice=$i" : ""
                        push!(seen, key)
                        P
                    end
                end
            end
        end
    end

    # --------------------------------------------------------
    # Plot observation partition stored in the quotient
    # --------------------------------------------------------
    if what == :obs_partition
        seen = Set{Int}()

        for (P, ob) in T.obs_partition
            if !isnothing(obs) && ob != obs
                continue
            end

            c = palette[mod1(ob + 2, length(palette))]

            lbl = if ob in seen || !show_labels
                ""
            elseif ob == -1
                "Terminal set"
            elseif ob == 0
                "Neutral region"
            else
                "Observation $ob"
            end
            push!(seen, ob)

            @series begin
                linecolor := c
                fillcolor := c
                fillalpha := fillalpha
                label := lbl
                P
            end
        end
    end
end

# @recipe function f(
#     T::PCBisimulationQuotient;
#     node = nothing,
#     slice = nothing,
#     obs = nothing,
#     by = :slice,
#     fillalpha = 0.25,
#     show_labels = false,
# )
#     palette = [:red, :blue, :green, :orange, :purple, :brown, :pink, :cyan]

#     for q in values(T.states)
#         if !isnothing(node) && q.node != node
#             continue
#         end
#         if !isnothing(slice) && q.slice != slice
#             continue
#         end
#         if !isnothing(obs) && q.obs != obs
#             continue
#         end

#         c = if by == :slice
#             palette[mod1(q.slice, length(palette))]
#         elseif by == :obs
#             palette[mod1(q.obs + 1, length(palette))]
#         elseif by == :node
#             palette[mod1(hash(q.node), length(palette))]
#         else
#             :blue
#         end

#         @series begin
#             linecolor := c
#             fillcolor := c
#             fillalpha := fillalpha
#             label := show_labels ? "q$(q.id)" : ""
#             q.set
#         end
#     end
# end
