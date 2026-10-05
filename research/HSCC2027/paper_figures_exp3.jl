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
# and adds the stacked pair, which the experiment does not draw.
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
include(joinpath(@__DIR__, "paper_hist.jl"))
include(joinpath(@__DIR__, "paper_stack.jl"))

has_latex() || error("needs pdflatex and pdftocairo on PATH")

figdir = joinpath(@__DIR__, "figures")
arms(r) = (r.a1.s.faces, r.a2.s.faces)

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
for r in (primal, dual)
    out = joinpath(figdir, "exp3_$(r.tag)_facets_histogram.png")
    hist_axes(facets_hist(arms(r)...), out)
    println("wrote ", basename(out))
end

# ── the two histograms stacked ───────────────────────────────────────────────────────────────────
# The pair above stand alone; this one is for reading the graphs against each other. Shared x and y
# limits and a single x axis, so a bar at a given height means the same in both rows, and the whole
# comparison fits one column where two separate figures would not.
let out = joinpath(figdir, "exp3_facets_histogram_stacked.png"),
    top = widest(arms(primal)..., arms(dual)...)

    latex_axes_stacked(
        [facets_hist(arms(primal)...; top), facets_hist(arms(dual)...; top)],
        out;
        rows = ["primal", "dual"],
        xlab = raw"\text{facets per cell}",
        ylab = raw"\text{cells}",
        panel_pt = 92,
        ylog = true,
        defs = HIST_DEFS,
        legend = HIST_LEGEND,
        # the dual's cells are far simpler, so its right half of the shared axis is empty and the
        # legend costs nothing there; over the primal it would sit on the worst-cell line
        legend_row = 2,
        # a column both panels leave empty: right of the tall first bars, clear of the dashed lines
        row_at = (0.24, 0.95),
    )
    println("wrote ", basename(out))
end

println("\nall five LaTeX-typeset figures written")
