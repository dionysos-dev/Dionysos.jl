# Where the geometric complexity goes.
#
# A path-complete certificate against the common Lyapunov function induced by it
# (`build_common_lyapunov`): the same certificate content, collapsed onto one node. The system,
# tolerances and observation regions are identical, so the only difference is whether the geometry
# sits in the template or in the memory.
#
# The total is essentially conserved — as the facet-budget bound requires, since it is bounded
# below by a quantity depending only on the system. What changes is the distribution, and that is
# what governs cost: the geometric operations of the refinement (pre-image, intersection, set
# difference) are superlinear in the facet count of the cells involved, so the same total spread
# over simpler cells is cheaper to manipulate.
#
# Produces the figure comparing the two distributions, written to JCVD/Figures/.

include(joinpath(@__DIR__, "common.jl"))

using Statistics
using Printf

const PCLF_FILE = joinpath(@__DIR__, "pclf_case.jld2")
const CLF_FILE = joinpath(@__DIR__, "clf_case.jld2")
const FIGURE_DIR = joinpath(@__DIR__, "JCVD", "Figures")

function load_quotient(filename)
    isfile(filename) || error(
        "$(basename(filename)) missing — run pclf_vs_clf.jl first to build and cache the quotients.",
    )
    optimizer = import_optimizer_jld2(filename)
    return MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
end

bisim_pclf = load_quotient(PCLF_FILE)
bisim_clf = load_quotient(CLF_FILE)

parts_pclf, facets_pclf = PCQ.cell_complexities(bisim_pclf)
parts_clf, facets_clf = PCQ.cell_complexities(bisim_clf)

const QUANTILES = [0.5, 0.9, 0.99, 1.0]

println(
    rpad("", 12),
    rpad("cells", 8),
    rpad("Σ", 10),
    rpad("mean", 9),
    join([rpad("q$(q)", 9) for q in QUANTILES]),
)
println("-"^76)
for (name, values) in [
    ("PCLF pieces", parts_pclf),
    ("CLF pieces", parts_clf),
    ("PCLF facets", facets_pclf),
    ("CLF facets", facets_clf),
]
    println(
        rpad(name, 12),
        rpad(length(values), 8),
        rpad(sum(values), 10),
        rpad(round(mean(values); digits = 2), 9),
        join([rpad(round(quantile(values, q); digits = 1), 9) for q in QUANTILES]),
    )
end

@printf(
    "\nTotal facets: PCLF %d, CLF %d — ratio %.3f\n",
    sum(facets_pclf),
    sum(facets_clf),
    sum(facets_clf) / sum(facets_pclf),
)
@printf(
    "Worst cell:   PCLF %d facets, CLF %d facets — ratio %.1f×\n",
    maximum(facets_pclf),
    maximum(facets_clf),
    maximum(facets_clf) / maximum(facets_pclf),
)

# ---------------------------------------------------------
# Figure
# ---------------------------------------------------------
# Log counts: most cells are simple and a long tail is not, and the tail is the point.

gr()

function distribution_panel(pclf_values, clf_values; bins, xlabel, title)
    fig = histogram(
        pclf_values;
        bins = bins,
        yscale = :log10,
        alpha = 0.55,
        label = "PCLF (4 nodes)",
        color = :steelblue,
        linecolor = :steelblue,
        xlabel = xlabel,
        ylabel = "number of cells",
        title = title,
        ylims = (0.7, 2e4),
        legend = :topright,
        framestyle = :box,
        grid = false,
    )
    histogram!(
        fig,
        clf_values;
        bins = bins,
        yscale = :log10,
        alpha = 0.55,
        label = "induced CLF (1 node)",
        color = :firebrick,
        linecolor = :firebrick,
    )
    # The two maxima, which is where the separation actually lives.
    vline!(
        fig,
        [maximum(pclf_values)];
        color = :steelblue,
        ls = :dash,
        lw = 1.5,
        label = "",
    )
    vline!(fig, [maximum(clf_values)]; color = :firebrick, ls = :dash, lw = 1.5, label = "")
    return fig
end

pieces_panel = distribution_panel(
    parts_pclf,
    parts_clf;
    bins = 0:4:450,
    xlabel = "polytopes per cell",
    title = "Pieces per semilinear cell",
)
facets_panel = distribution_panel(
    facets_pclf,
    facets_clf;
    bins = 0:20:1800,
    xlabel = "facets per cell",
    title = "Facets per semilinear cell",
)

fig = plot(
    pieces_panel,
    facets_panel;
    layout = (1, 2),
    size = (900, 340),
    left_margin = 5Plots.mm,
    bottom_margin = 5Plots.mm,
    top_margin = 3Plots.mm,
)

mkpath(FIGURE_DIR)
savefig(fig, joinpath(FIGURE_DIR, "complexity_distribution.pdf"))
savefig(fig, joinpath(FIGURE_DIR, "complexity_distribution.png"))
println("\nwrote $(joinpath("JCVD", "Figures", "complexity_distribution")).{pdf,png}")
display(fig)
