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

const PALETTE = [
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

# Accumulate the polygons of `parts` into one coordinate pair, separated by NaN.
#
# A quotient has thousands of cells, and Plots' cost per series grows with the series already in
# the figure, so one series per cell is quadratic: measured 0.31 ms per cell at 250 cells against
# 4.25 ms at 8000. NaN breaks both the fill and the outline, so a merged series draws each cell
# individually — a single-colour figure comes out byte-identical to the per-cell version.
function _append_polygons!(X::Vector{Float64}, Y::Vector{Float64}, parts)
    for P in parts
        x, y = LazySets.plot_recipe(P)
        isempty(x) && continue
        append!(X, x)
        append!(Y, y)
        push!(X, NaN)
        push!(Y, NaN)
    end
    return nothing
end

# Bucket `(colour, parts)` pairs into one coordinate pair per colour, keeping first-seen order so
# the drawing order of the colours is the order the cells were visited in.
function _group_by_colour(coloured_parts)
    groups = Tuple{Any, Vector{Float64}, Vector{Float64}}[]
    index = Dict{Any, Int}()
    for (c, parts) in coloured_parts
        gi = get!(index, c) do
            push!(groups, (c, Float64[], Float64[]))
            return length(groups)
        end
        _append_polygons!(groups[gi][2], groups[gi][3], parts)
    end
    return groups
end

# ------------------------------------------------------------
# Full quotient recipe
# ------------------------------------------------------------
# `merge_series` trades per-cell series for per-colour series. It is off automatically when
# `show_labels` is on, because a legend entry per abstract state needs a series per abstract
# state. Merging also fixes the drawing order to colour-group order rather than cell order,
# which is invisible for disjoint cells but does decide which of two *overlapping* cells is on
# top — pass `merge_series = false` when that matters.
@recipe function f(
    quotient::PCBisimulationQuotient;
    what = :states,
    node = nothing,
    slice = nothing,
    obs = nothing,
    state_ids = nothing,
    by = :state,
    fillalpha = 0.25,
    linewidth = 1.5,
    seriesalpha = 0.9,
    show_labels = false,
    show_contours = false,
    user_color = nothing,
    merge_series = true,
)
    palette = PALETTE

    local_linealpha = show_contours ? 1.0 : 0.0
    state_id_set = isnothing(state_ids) ? nothing : Set(state_ids)
    merged = merge_series && !show_labels

    if what == :states
        qlist = [
            q for q in values(quotient.states) if (isnothing(node) || q.node == node) &&
                (isnothing(slice) || q.slice == slice) &&
                (isnothing(obs) || q.obs == obs) &&
                (isnothing(state_id_set) || q.id in state_id_set)
        ]

        sort!(qlist; by = q -> q.id)

        colour_of(k, q) =
            if !isnothing(user_color)
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

        if merged
            for (c, X, Y) in _group_by_colour(
                (colour_of(k, q), q.set.array) for (k, q) in enumerate(qlist)
            )
                isempty(X) && continue
                @series begin
                    seriestype := :shape
                    fillcolor := c
                    linecolor := c
                    fillalpha := fillalpha
                    linealpha := local_linealpha
                    linewidth := linewidth
                    seriesalpha := seriesalpha
                    label := ""
                    X, Y
                end
            end
        else
            for (k, q) in enumerate(qlist)
                c = colour_of(k, q)
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
    end

    if what == :slices
        groups = Tuple{Any, Int, UT.SemiLinearSet}[]
        for (nd, slice_list) in quotient.slices
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

        if merged
            for (c, X, Y) in _group_by_colour(
                (palette[mod1(i, length(palette))], S.array) for (_, i, S) in groups
            )
                isempty(X) && continue
                @series begin
                    seriestype := :shape
                    fillcolor := c
                    linecolor := c
                    fillalpha := fillalpha
                    linealpha := local_linealpha
                    linewidth := linewidth
                    seriesalpha := seriesalpha
                    label := ""
                    X, Y
                end
            end
        else
            seen = Set{Tuple{Any, Int}}()
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
                            (show_labels && !(key in seen) && j == 1) ?
                            "node=$nd, slice=$i" : ""
                        P
                    end
                end

                push!(seen, key)
            end
        end
    end
end
