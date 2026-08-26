# ============================================================
# Plot recipes for the bisimulation quotient
# ============================================================

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
