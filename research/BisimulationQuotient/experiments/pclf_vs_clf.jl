# Path-complete versus common Lyapunov function: how much cheaper are the quotient's cells?
#
# Both quotients are built for the same problem. One uses the path-complete family over a
# four-node observer graph; the other uses the common Lyapunov function induced by it, so the
# certificates are directly comparable. The figure histograms the per-cell complexity — pieces
# and faces — because that is what decides whether the quotient is workable.
#
# Building both takes minutes, so each is cached in a `*.jld2` beside this script and reused.

include(joinpath(dirname(@__DIR__), "common.jl"))

gr()

const PCLF_FILE = joinpath(@__DIR__, "pclf_case.jld2")
const CLF_FILE = joinpath(@__DIR__, "clf_case.jld2")

(; f, problem) = observer_graph_problem(; p = 1.7)

"""
Build both quotients and cache them, unless the caches already exist.

Pass `force = true` after changing the problem or the tolerances, or the stale caches are
loaded instead.
"""
function compute_and_cache(;
    force::Bool = false,
    atol = 1e-4,
    level_tol = 1e-2,
    max_slices = 20,
)
    if !force && isfile(PCLF_FILE) && isfile(CLF_FILE)
        println("Reusing cached quotients; pass force = true to rebuild.")
        return nothing
    end

    pclf_poly = observer_graph_pclf(f)
    println("Computed JSR upper bound / contraction rate = ", pclf_poly.JSRapprox)

    println("--- path-complete certificate ---")
    (; optimizer) = build_quotient(
        problem,
        pclf_poly;
        atol = atol,
        level_tol = level_tol,
        max_slices = max_slices,
    )
    export_optimizer_jld2(optimizer, PCLF_FILE)

    println("--- induced common certificate ---")
    clf_poly = PCLF.build_common_lyapunov(pclf_poly)
    (; optimizer) = build_quotient(
        problem,
        clf_poly;
        atol = atol,
        level_tol = level_tol,
        max_slices = max_slices,
    )
    export_optimizer_jld2(optimizer, CLF_FILE)
    return nothing
end

function load_quotient(filename)
    isfile(filename) || error("$filename missing -- run compute_and_cache() first.")
    optimizer = import_optimizer_jld2(filename)
    return MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
end

compute_and_cache()

bisim_pclf = load_quotient(PCLF_FILE)
bisim_clf = load_quotient(CLF_FILE)

n_parts_pclf, n_faces_pclf = PCQ.cell_complexities(bisim_pclf)
n_parts_clf, n_faces_clf = PCQ.cell_complexities(bisim_clf)

# Log counts: most cells are simple and a long tail is not, which is the whole point.
const LOG_TICKS = ([1, 10, 100, 1000], ["10⁰", "10¹", "10²", "10³"])

function complexity_histogram(pclf_values, clf_values; bins, xlabel, title, xlims)
    fig = histogram(
        pclf_values;
        bins = bins,
        xlabel = xlabel,
        ylabel = "# cells",
        title = title,
        yscale = :log10,
        xlims = xlims,
        ylims = (:auto, 1e4),
        yticks = LOG_TICKS,
        legend = :topright,
        alpha = 0.45,
        label = "PCLF",
    )
    histogram!(fig, clf_values; bins = bins, alpha = 0.45, label = "CLF", yscale = :log10)
    return fig
end

p1 = complexity_histogram(
    n_parts_pclf,
    n_parts_clf;
    bins = 0:2:450,
    xlabel = "# polytopes",
    title = "Pieces per semilinear cell",
    xlims = (0, 450),
)
p2 = complexity_histogram(
    n_faces_pclf,
    n_faces_clf;
    bins = 0:10:2000,
    xlabel = "# faces",
    title = "Faces per semilinear cell",
    xlims = (0, 2000),
)

fig = plot(
    p1,
    p2;
    layout = (1, 2),
    size = (1400, 500),
    margin = 8Plots.mm,
    titlefontsize = 14,
    guidefontsize = 12,
    tickfontsize = 10,
)
display(fig)
