# Every experiment 3 figure that LaTeX typesets, in the paper's own fonts.
#
#     julia --project=. paper_figures_exp3.jl          [LEVELS=n] [CACHE=0]
#
# Run after `experiment3_diverse_pieces.jl`, which draws the same figures through GR alone — fast,
# self-contained, no LaTeX, and therefore in Computer Modern rather than acmart's Libertine. This
# script replaces four of the six:
#
#     exp3_{primal,dual}_approach1_quotient      exp3_{primal,dual}_facets_histogram
#
# The two 3-D figures are CairoMakie and cannot take this route; the experiment gives them Libertine
# directly, through Makie's own font loading.
#
# It runs with CACHE=1 and TIMING=0, so the quotients come off disk and the four timed builds per
# graph are skipped: about a minute here against ten for a full experiment run. Nothing it prints is
# a measurement.
#
# Needs `pdflatex` and `pdftocairo` on PATH. Rebuilding the experiment's own figures needs neither.

ENV["FIGURES"] = "0"
ENV["TIMING"] = "0"           # the timed builds are the measurement, and no figure needs them
ENV["CACHE"] = get(ENV, "CACHE", "1")
include(joinpath(@__DIR__, "experiment3_diverse_pieces.jl"))
include(joinpath(@__DIR__, "paper_fonts.jl"))
include(joinpath(@__DIR__, "paper_stack.jl"))

has_latex() || error("needs pdflatex and pdftocairo on PATH")

"Widest cell of either arm, over the graphs given — the binning two panels have to share."
widest(rs...) = maximum(max(maximum(r.a1.s.faces), maximum(r.a2.s.faces)) for r in rs)

"""
The facets-per-cell histogram for one graph, bars only — the axes come from LaTeX.

`top` is where the bins stop. It defaults to this graph's own widest cell, which is what a figure
standing on its own wants; panels meant to be read against each other pass the pair's widest instead,
so that a bar of a given width means the same in both.
"""
function facets_hist(r; top = widest(r))
    bins = range(0, 1.02 * top; length = 45)
    h = histogram(
        r.a1.s.faces;
        bins,
        yscale = :log10,
        fillrange = 0.7,
        alpha = 0.55,
        color = :firebrick,
        linecolor = :firebrick,
        label = "",
        legend = false,
        ylims = (0.7, 5e3),
    )
    histogram!(
        h,
        r.a2.s.faces;
        bins,
        yscale = :log10,
        fillrange = 0.7,
        alpha = 0.55,
        color = :steelblue,
        linecolor = :steelblue,
        label = "",
    )
    vline!(h, [r.a1.s.max_f]; color = :firebrick, ls = :dash, lw = 1.8, label = "")
    vline!(h, [r.a2.s.max_f]; color = :steelblue, ls = :dash, lw = 1.8, label = "")
    return h
end

# Julia's :firebrick and :steelblue, as RGB, so the legend swatches match the bars
const HIST_DEFS = [
    raw"\definecolor{firebrick}{RGB}{178,34,34}",
    raw"\definecolor{steelblue}{RGB}{70,130,180}",
]
const HIST_LEGEND =
    [("firebrick, line width=2pt", "iCLF"), ("steelblue, line width=2pt", "PCLF")]
figdir = joinpath(@__DIR__, "figures")

# ── approach 1's quotient, one per graph ─────────────────────────────────────────────────────────
for r in (primal, dual)
    p = plot(; aspect_ratio = :equal, legend = false)
    plot!(
        p,
        r.a1.b.quotient;
        what = :states,
        by = :state,
        show_contours = true,
        linewidth = 0.2,
        fillalpha = 0.9,
    )
    out = joinpath(figdir, "exp3_$(r.tag)_approach1_quotient.png")
    latex_axes(p, out)
    println("wrote ", basename(out))
end

# ── facets per cell, one per graph ───────────────────────────────────────────────────────────────
# Most cells are simple and the construction pays for the tail, so the count is logarithmic.
# `fillrange` pins both series to one baseline, which a log axis otherwise picks per series.
for r in (primal, dual)
    out = joinpath(figdir, "exp3_$(r.tag)_facets_histogram.png")
    latex_axes(
        facets_hist(r),
        out;
        xlab = raw"\text{facets per cell}",
        ylab = raw"\text{cells}",
        height_pt = 105,
        # a fifth the height of the square panels, so the shared 11/14 pt is far too big here;
        # the legend takes the tick size too
        tick_pt = 9,
        label_pt = 10,
        ylog = true,
        defs = HIST_DEFS,
        # The dashed lines mark each arm's worst cell, 510 against 64 on the primal graph. They are
        # deliberately absent from the legend, which keeps it to two entries; the caption carries
        # the numbers, and the script prints them.
        legend = HIST_LEGEND,
    )
    println("wrote ", basename(out))
end

# ── the two histograms stacked ───────────────────────────────────────────────────────────────────
# The pair above stand alone; this one is for reading the graphs against each other. Shared x and y
# limits and a single x axis, so a bar at a given height means the same in both rows, and the whole
# comparison fits one column where two separate figures would not.
let out = joinpath(figdir, "exp3_facets_histogram_stacked.png"), top = widest(primal, dual)
    latex_axes_stacked(
        [facets_hist(primal; top), facets_hist(dual; top)],
        out;
        rows = ["primal", "dual"],
        xlab = raw"\text{facets per cell}",
        ylab = raw"\text{cells}",
        panel_pt = 92,
        ylog = true,
        defs = HIST_DEFS,
        legend = HIST_LEGEND,
        # the dual reaches half as far along the shared x axis, so its right half is empty and the
        # legend costs nothing there; over the primal it would sit on the 510-facet line
        legend_row = 2,
        # a column both panels leave empty: right of the tall first bars, left of the dual's
        # dashed lines, and nowhere near the legend in the row below
        row_at = (0.24, 0.95),
    )
    println("wrote ", basename(out))
end

println("\nall five LaTeX-typeset figures written")
