# Every experiment 1 figure that LaTeX typesets, in the paper's own fonts.
#
#     julia --project=. paper_figures/paper_figures_exp1.jl   [CACHE=0]
#
# Run after `experiment1_gol_lazar_belta.jl`, which draws the same figures through GR alone — fast,
# self-contained, no LaTeX, and therefore in Computer Modern rather than acmart's Libertine. This
# script replaces seven of them:
#
#     exp1_problem              exp1_approach1_quotient      exp1_facets_histogram
#     exp1_approach{1,2}_synthesis     exp1_approach{1,2}_verification
#
# The two 3-D figures are CairoMakie and cannot take this route; the experiment gives them Libertine
# directly, through Makie's own font loading.
#
# Needs `pdflatex` and `pdftocairo` on PATH. Rebuilding the experiment's own figures needs neither.
#
# WHY ONE SCRIPT AND NOT FIVE. Loading the experiment costs about three minutes — packages, the
# certificate, 35 MB of cached quotients, both co-safe solves and the witness tree — against roughly
# a second of LaTeX per figure. Doing it once instead of once per figure is the whole difference
# between four minutes and twelve, and it collapses the ordering rule to "experiment, then this".

ENV["FIGURES"] = "0"          # the witness runs are computed outside that guard
ENV["TIMING"] = "0"           # the timed builds are the measurement, and no figure needs them
ENV["CACHE"] = get(ENV, "CACHE", "1")
include(joinpath(dirname(@__DIR__), "experiment1_gol_lazar_belta.jl"))
include(joinpath(@__DIR__, "paper_fonts.jl"))
include(joinpath(@__DIR__, "paper_hist.jl"))

has_latex() || error("needs pdflatex and pdftocairo on PATH")
d = joinpath(dirname(@__DIR__), "figures")

# ── the sets the experiment's skipped figure block would have built ──────────────────────────────
τD = ΓX * γ^(NB - 1)
D_pieces = [UT._as_hpolytope(PCLF.get_sublevel_set(pclf.pieces[nd], τD)) for nd in nodes]
D1 = reduce(LazySets.intersection, D_pieces)
X_pieces = [UT._as_hpolytope(PCLF.get_sublevel_set(pclf.pieces[nd], ΓX)) for nd in nodes]
X1 = reduce(LazySets.intersection, X_pieces)

const PROBLEM_COLOURS = [:darkgreen, :navy, :darkorange]
const LBL = 16
# R_2 is too narrow to hold a 16 pt label, so its label sits beside the region rather than in it.
LABEL_NUDGE = Dict("R1" => [0.9, 0.0], "R2" => [2.6, -1.6])

# ══ the GR-rastered figures ══════════════════════════════════════════════════════════════════════
# These carry thousands of polytopes, far past what pgfplots will compile, so GR draws the graphics
# and LaTeX draws only the frame, the ticks, the labels and the legend. `gr()` is still the active
# backend here, from the experiment.

# ── approach 1's quotient ────────────────────────────────────────────────────────────────────────
let p = plot(; aspect_ratio = :equal, legend = false)
    plot!(
        p,
        a1.built.quotient;
        what = :states,
        by = :state,
        show_contours = true,
        linewidth = 0.2,
        fillalpha = 0.9,
    )
    # `show_contours = true` leaves a subplot-level line alpha of 0.5 that later series inherit,
    # which is what makes these outlines vanish if it is not pinned back.
    for (i, (_, R)) in enumerate(REGIONS)
        plot!(
            p,
            R;
            fillalpha = 0.0,
            linecolor = REGION_COLOURS[i],
            linealpha = 1.0,
            linewidth = 1.4,
            label = "",
        )
    end
    out = joinpath(d, "exp1_approach1_quotient.png")
    latex_axes(p, out; height_pt = 272)          # the same axis box as exp1_problem
    println("wrote ", basename(out))
end

# ── the four specification panels ────────────────────────────────────────────────────────────────
"The certified set of one arm under one quantifier, with the runs that witness it."
function spec_panel(arm, res, runs)
    p = plot(; aspect_ratio = :equal, legend = false)
    for (ids, c) in ((res.lost, LOST), (res.won, WON))
        isempty(ids) && continue
        plot!(
            p,
            arm.built.quotient;
            what = :states,
            state_ids = collect(ids),
            show_contours = false,
            user_color = c,
            fillalpha = 1.0,
        )
    end
    # the regions go back on top: the runs cross them, and a boundary the eye cannot follow stops
    # being a reference
    for (i, (_, R)) in enumerate(REGIONS)
        plot!(
            p,
            R;
            fillalpha = 0.0,
            linecolor = REGION_COLOURS[i],
            linealpha = 1.0,
            linewidth = 1.4,
            label = "",
        )
    end
    for P in (arm.built.D isa LazySets.UnionSetArray ? arm.built.D.array : [arm.built.D])
        plot!(
            p,
            P;
            fillalpha = 0.0,
            linecolor = :black,
            linestyle = :dash,
            linealpha = 1.0,
            linewidth = 0.8,
            label = "",
        )
    end
    for (X, col) in runs
        plot!(
            p,
            first.(X),
            last.(X);
            color = col,
            linealpha = 1.0,
            linewidth = 1.2,
            marker = :circle,
            markersize = 2.0,
            markerstrokewidth = 0,
            label = "",
        )
    end
    scatter!(
        p,
        [xw[1]],
        [xw[2]];
        color = :white,
        markerstrokecolor = :black,
        markerstrokewidth = 0.8,
        markersize = 3,
        label = "",
    )
    return p
end

# pgfplots styles matching the run colours, for the LaTeX-typeset legend
sty(
    col,
) = "$col, line width=1.2pt, mark=*, mark size=1.1pt, mark options={solid, fill=$col}"

for (name, arm, sim) in (("exp1_approach1", a1, sim1), ("exp1_approach2", a2, sim2))
    out = joinpath(d, "$(name)_synthesis.png")
    latex_axes(
        spec_panel(arm, arm.syn, [(sim.X, RUN_COLOURS[1])]),
        out;
        height_pt = 272,
        legend = [(sty("black"), "\$\\sigma = $(word_label(Int.(sim.U)))\$")],
    )
    println("wrote ", basename(out))

    out = joinpath(d, "$(name)_verification.png")
    latex_axes(
        spec_panel(arm, arm.ver, [(b.X, RUN_COLOURS[i]) for (i, b) in enumerate(shown)]),
        out;
        height_pt = 272,
        legend = [
            (sty(string(RUN_COLOURS[i])), "\$\\sigma_$i = $(word_label(b.word))\$") for
            (i, b) in enumerate(shown)
        ],
    )
    println("wrote ", basename(out))
end

# ── facets per cell ──────────────────────────────────────────────────────────────────────────────
# This quotient's tail runs an order of magnitude higher than experiments 2 and 3, so the count axis
# goes to 5e4 rather than their 5e3.
let out = joinpath(d, "exp1_facets_histogram.png")
    hist_axes(facets_hist(a1.s.faces, a2.s.faces; ymax = 5e4), out)
    println("wrote ", basename(out))
end

# ══ the problem figure, drawn wholly by LaTeX ════════════════════════════════════════════════════
# A dozen polygons, so pgfplots can compile the geometry itself — which also puts the `\mathcal`
# labels in newtxmath rather than in the raster. Backend switched last, since everything above
# needs GR.
using PGFPlotsX
pgfplotsx()

# Exactly what acmart loads under pdflatex; see its font block. Under lualatex or xelatex the class
# takes `\mathcal` from Latin Modern instead, which would NOT match a pdflatex-built paper.
PGFPlotsX.CUSTOM_PREAMBLE = [
    raw"\usepackage[T1]{fontenc}",
    raw"\usepackage[tt=false, type1=true]{libertine}",
    raw"\usepackage[libertine]{newtxmath}",
]

"Vertices in cyclic order, first repeated: `vertices_list` promises no order, and pgfplots strokes
the path as given without closing it, so an open ring loses its last edge."
function ring(P)
    V = LazySets.vertices_list(UT._as_hpolytope(P))
    c = sum(V) / length(V)
    V = V[sortperm([atan(v[2] - c[2], v[1] - c[1]) for v in V])]
    push!(V, V[1])
    return Shape(first.(V), last.(V))
end

let p = plot(;
        aspect_ratio = :equal,
        legend = false,
        xlabel = L"x_1",
        ylabel = L"x_2",
        tickfontsize = 11,
        guidefontsize = 14,
        grid = false,
        framestyle = :axes,
        size = (241, 272),
        # pgfplots sizes the axis box itself and ignores Plots' `size`; width alone would stretch
        # the plot and destroy equal aspect, so the height carries the data's own ratio.
        #
        # The thin frame closes the axis box at the data limits, as the raster's frame does in the
        # `latex_axes` figures this one sits beside; the ticks stay on the left and bottom.
        extra_kwargs = Dict(
            :subplot => Dict(
                "width" => "241pt",
                "height" => "272pt",
                "scale only axis" => "true",
                "after end axis/.code" =>
                    raw"{\draw[line width=0.4pt] (rel axis cs:0,0) rectangle (rel axis cs:1,1);}",
            ),
        ),
    )
    for P in X_pieces
        plot!(
            p,
            ring(P);
            fillcolor = :grey93,
            fillalpha = 1.0,
            linecolor = :grey70,
            linewidth = 0.6,
            label = "",
        )
    end
    plot!(p, ring(X1); fillalpha = 0.0, linecolor = :black, linewidth = 1.1, label = "")
    for (i, (_, R)) in enumerate(REGIONS)
        plot!(
            p,
            ring(R);
            fillcolor = PROBLEM_COLOURS[i],
            fillalpha = 0.30,
            linecolor = PROBLEM_COLOURS[i],
            linewidth = 1.4,
            label = "",
        )
    end
    for P in D_pieces
        plot!(
            p,
            ring(P);
            fillcolor = :red,
            fillalpha = 0.40,
            linecolor = :red,
            linewidth = 0.7,
            label = "",
        )
    end
    plot!(
        p,
        ring(D1);
        fillcolor = :grey78,
        fillalpha = 1.0,
        linecolor = :black,
        linewidth = 0.8,
        label = "",
    )

    for (i, (nm, R)) in enumerate(REGIONS)
        c = cell_point(R) + get(LABEL_NUDGE, nm, [0.0, 0.0])
        annotate!(
            p,
            c[1],
            c[2],
            Plots.text(latexstring("R_", i), LBL, PROBLEM_COLOURS[i], :center),
        )
    end
    let c = cell_point(D1)
        annotate!(p, c[1], c[2], Plots.text(L"\mathcal{D}_1", LBL, :black, :center))
    end
    # placed by hand: taken from the outermost vertex it lands on its own red boundary and level
    # with R_3, as LABEL_NUDGE already does for R_1 and R_2
    annotate!(p, 4.1, 4.2, Plots.text(L"\mathcal{D}_2", LBL, :red, :center))
    let V = LazySets.vertices_list(X1)
        v = V[argmin([w[2] for w in V])]
        annotate!(
            p,
            0.88 * v[1],
            0.88 * v[2],
            Plots.text(L"\mathcal{X}_1", LBL, :black, :center),
        )
    end
    let V = reduce(vcat, LazySets.vertices_list.(X_pieces)),
        cs = LazySets.constraints_list(X1)

        far = filter(w -> maximum(dot(c.a, w) - c.b for c in cs) > 1e-8, V)
        cand = isempty(far) ? V : far
        v = cand[argmax([w[1] + w[2] for w in cand])]
        annotate!(
            p,
            1.07 * v[1] + 0.5,
            1.07 * v[2],
            Plots.text(L"\mathcal{X}_2", LBL, :black, :center),
        )
    end

    pdf = joinpath(d, "exp1_problem_paper.pdf")
    png = joinpath(d, "exp1_problem.png")
    savefig(p, pdf)
    Base.run(`pdftocairo -png -r 400 -singlefile $pdf $(splitext(png)[1])`)
    rm(pdf; force = true)
    crop_png!(png)
    println("wrote ", basename(png))
end

gr()                                   # leave the session on the backend it arrived with
println("\nall seven LaTeX-typeset figures written")
