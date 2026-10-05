# The facets-per-cell histogram, which all three experiments draw the same way.
#
# One cell of a quotient is a union of polytopes, and every geometric primitive of the construction
# is charged per facet, so this distribution is the cost argument itself. Both arms go on one pair of
# axes with a logarithmic count, since most cells are simple and the construction pays for the tail.
#
# Include after paper_fonts.jl: `hist_axes` is `latex_axes` with the settings they share.

# Julia's :firebrick and :steelblue, as RGB, so the legend swatches match the bars
const HIST_DEFS = [
    raw"\definecolor{firebrick}{RGB}{178,34,34}",
    raw"\definecolor{steelblue}{RGB}{70,130,180}",
]
# The dashed lines mark each arm's worst cell. They are deliberately absent, which keeps the legend
# to two entries; the caption carries the numbers, and the experiments print them.
const HIST_LEGEND =
    [("firebrick, line width=2pt", "iCLF"), ("steelblue, line width=2pt", "PCLF")]

"Widest cell among the face counts given — the binning panels read against each other must share."
widest(vs...) = maximum(maximum(v) for v in vs)

"""
    facets_hist(f1, f2; top, ymax)

The two arms' facet counts as one overlaid histogram, bars only — the axes come from `hist_axes`.
`fillrange` pins both series to a single baseline, which a log axis otherwise picks per series.

`top` is where the bins stop. It defaults to the widest cell of this pair, which is what a figure
standing on its own wants; panels meant to be compared pass the widest over all of them instead, so
that a bar of a given width means the same in every panel.
"""
function facets_hist(f1, f2; top = widest(f1, f2), ymax = 5e3)
    bins = range(0, 1.02 * top; length = 45)
    common =
        (; bins, yscale = :log10, fillrange = 0.7, alpha = 0.55, label = "", legend = false)
    h = histogram(
        f1;
        common...,
        color = :firebrick,
        linecolor = :firebrick,
        ylims = (0.7, ymax),
    )
    histogram!(h, f2; common..., color = :steelblue, linecolor = :steelblue)
    vline!(h, [maximum(f1)]; color = :firebrick, ls = :dash, lw = 1.8, label = "")
    vline!(h, [maximum(f2)]; color = :steelblue, ls = :dash, lw = 1.8, label = "")
    return h
end

"""
    hist_axes(h, out; kwargs...)

`latex_axes` with what every facets histogram wants: a fifth the height of the square panels, so the
shared 11/14 pt is far too big here, and the legend takes the tick size too.
"""
function hist_axes(h, out; kwargs...)
    return latex_axes(
        h,
        out;
        xlab = raw"\text{facets per cell}",
        ylab = raw"\text{cells}",
        height_pt = 105,
        tick_pt = 9,
        label_pt = 10,
        ylog = true,
        defs = HIST_DEFS,
        legend = HIST_LEGEND,
        kwargs...,
    )
end
