# EXPERIMENT 2 — the same system on the primal and on the dual De Bruijn graph.
#
#     julia --project=. experiment2_primal_and_dual.jl   [SAMPLES=n] [FIGURES=0] [CACHE=0] [TIMING=0]
#
# Self-contained: this file defines its own system, certificate, levels and figures.
#
# One system, one certificate family, one working set. The only thing that changes is what a graph
# node means:
#
#   primal   the node records the mode just PLAYED      — every node emits every mode (complete)
#   dual     the node commits to the mode played NEXT   — every node emits one mode (co-complete)
#
# Approach 2 lifts the abstraction over the graph's nodes, so the obvious expectation is that it
# always builds the bigger quotient. It does not: on the primal graph it builds MORE cells than
# approach 1, on the dual FEWER, and it is faster on both. The cell count does not predict the cost;
# what governs the reversal is left to future work.
#
# The node pieces are kept close here, by sharing one template between the nodes, which is what
# isolates the graph. Pieces that separate help approach 2 on the primal graph and hurt it on the
# dual — that is experiment 3.

using StaticArrays, JuMP, Clarabel, LazySets, Plots, Printf, LaTeXStrings
import HybridSystems, Statistics, PNGFiles, JLD2
import Plots.PlotMeasures: mm
import MathOptInterface as MOI

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const PCQ = DI.Optim.Abstraction.PCLFBisimulationQuotient
const PCLF = UT.PathCompleteFramework

gr()

const SAMPLES = parse(Int, get(ENV, "SAMPLES", "1"))
const FIGURES = get(ENV, "FIGURES", "1") == "1"
# Each graph builds six quotients: one per arm, and a warm-up plus SAMPLES timed rounds per arm.
# CACHE=1 keeps the first two on disk, which is what makes redrawing the figures cheap; TIMING=0
# drops the other four, which are the measurement and are never cached.
const USE_CACHE = get(ENV, "CACHE", "1") == "1"
const TIMING = get(ENV, "TIMING", "1") == "1"
const ATOL = 1e-3
const THETA = π / 6        # the shared template's orientation
const SOLVER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)
const REGION_COLOURS = [:black, :navy]

# ── the problem ──────────────────────────────────────────────────────────────────────────────────
A1 = @SMatrix [0.70 0.10; 0.00 0.65]
A2 = @SMatrix [0.60 -0.15; 0.10 0.55]
# The modes are declared the CONTROLLER's. `discreteswitchedsystem` would otherwise make them
# autonomous, i.e. the environment's.
f = ST.with_switching(
    HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)]),
    HybridSystems.ControlledSwitching(),
)
X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
R1 = LazySets.Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
R2 = LazySets.Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])
const REGIONS = [("R1", R1), ("R2", R2)]
const REGIONS_H = [UT._as_hpolytope(R) for (_, R) in REGIONS]
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])

println("="^84, "\nEXPERIMENT 2 — one system, the primal and the dual graph\n", "="^84)

# ── the two approaches, on one graph ─────────────────────────────────────────────────────────────
function build(cert, nb_levels, ΓX)
    opt = MOI.instantiate(PCQ.OptimizerBisimulationQuotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("pclf"), cert)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(opt, MOI.RawOptimizerAttribute("atol"), ATOL)
    MOI.set(opt, MOI.RawOptimizerAttribute("nb_levels"), nb_levels)
    MOI.set(opt, MOI.RawOptimizerAttribute("max_slices"), nb_levels)
    MOI.set(opt, MOI.RawOptimizerAttribute("ΓX"), ΓX)
    MOI.optimize!(opt)
    return (;
        quotient = MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
        D = MOI.get(opt, MOI.RawOptimizerAttribute("D")),
    )
end

"What `cert` builds, from cache/ when it is there. Timed builds never come through here."
function build_cached(name, tag, cert, nb_levels, ΓX)
    path = joinpath(@__DIR__, "cache", "exp2_$(tag)_$(name)_$(nb_levels).jld2")
    if USE_CACHE && isfile(path)
        return JLD2.jldopen(path, "r") do file
            return (; quotient = file["quotient"], D = file["D"])
        end
    end
    b = build(cert, nb_levels, ΓX)
    if USE_CACHE
        mkpath(dirname(path))
        JLD2.jldopen(path, "w") do file
            file["quotient"] = b.quotient
            return file["D"] = b.D
        end
    end
    return b
end

function timed(g, n)
    g()
    ts = Float64[]
    for _ in 1:n
        GC.gc()
        t0 = time_ns()
        g()
        push!(ts, (time_ns() - t0) / 1e9)
    end
    return Statistics.median(ts)
end

function stats(q)
    parts, faces = PCQ.cell_complexities(q)
    return (;
        cells = length(parts),
        parts,
        faces,
        max_p = maximum(parts),
        max_f = maximum(faces),
        mean_p = Statistics.mean(parts),
        mean_f = Statistics.mean(faces),
    )
end

function run_graph(dual)
    label = dual ? "DUAL (co-complete)" : "PRIMAL (complete)"
    tag = dual ? "dual" : "primal"
    graph = PCLF.generate_DeBruijn_edges(2, 1; dual = dual)
    nodes = sort(collect(graph.verts); by = string)
    Rot = [cos(THETA) -sin(THETA); sin(THETA) cos(THETA)]
    pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
        f,
        graph,
        SOLVER;
        Gmats = Dict(nd => Rot for nd in nodes),
        MLF = true,
        verbose = false,
    )
    isinf(pclf.JSRapprox) && error("no certificate on the $label graph")
    common = PCLF.build_common_lyapunov(pclf)
    γ = pclf.JSRapprox

    # The level at which the induced COMMON covers the working set, imposed on both approaches. The
    # certificate keeps its natural scale: rescaling would change the row norms, and with them the
    # geometric thickness of `atol`.
    ΓX = PCQ.compute_tau_X(common, UT._as_hpolytope(X))
    # The construction's own stopping rule: the first level whose terminal set no longer meets a
    # region, for both certificates.
    clears(j) = all(
        c -> PCQ.all_nodes_clear_regions(c, ΓX * γ^(j - 1), REGIONS_H; tol = 1e-2),
        (pclf, common),
    )
    NB = something(findfirst(clears, 1:60), -1)
    NB > 0 || error("no region-free terminal level on the $label graph")

    common_parts = let S = PCLF.get_sublevel_set(common.pieces[:clf], 1.0; atol = 1e-6)
        S isa LazySets.UnionSetArray ? length(S.array) : 1
    end

    println("\n", "#"^84, "\n$label De Bruijn, order 1\n", "#"^84)
    @printf(
        "complete: %-5s  co-complete: %-5s  certified rate %.6f\n",
        PCLF.is_complete(graph, 1:2),
        PCLF.is_co_complete(graph, 1:2),
        γ
    )
    @printf(
        "induced common: %d convex part(s)   levels: ΓX = %.4f, %d of them\n",
        common_parts,
        ΓX,
        NB
    )
    # A polyhedral piece is the gauge `V_s(x) = max_i |(G x)_i| / w_i`, whose Γ-sublevel set is the
    # symmetric polytope `{x : |G x| ≤ Γ w}`.
    println("  the certificate, one polyhedral piece per node:")
    for nd in nodes
        pc = pclf.pieces[nd]
        @printf("    node %-6s\n", string(nd))
        for i in 1:size(pc.G, 1)
            @printf(
                "      G[%d,:] = %9.4f %9.4f      w[%d] = %8.4f\n",
                i,
                pc.G[i, 1],
                pc.G[i, 2],
                i,
                pc.w[i]
            )
        end
    end

    # Build both arms before timing either. A timed region pays garbage collection proportional to
    # the live heap, so timing an arm while the other does not yet exist favours whichever goes
    # first, by enough to reverse the verdict.
    function arm(name, cert)
        b = build_cached(name, tag, cert, NB, ΓX)
        parts = b.D isa LazySets.UnionSetArray ? b.D.array : [UT._as_hpolytope(b.D)]
        for (nm, R) in REGIONS, P in parts
            LazySets.isdisjoint(P, UT._as_hpolytope(R)) ||
                error("[$label / $name] the terminal set D meets region $nm")
        end
        return (; b, s = stats(b.quotient), cert)
    end
    a1, a2 = arm("approach1", common), arm("approach2", pclf)
    time_arm(a) = TIMING ? timed(() -> build(a.cert, NB, ΓX), SAMPLES) : NaN
    a1, a2 = (; a1..., t = time_arm(a1)), (; a2..., t = time_arm(a2))
    for (name, a) in (("approach 1", a1), ("approach 2", a2))
        @printf(
            "  %-10s %10s  %5d cells  %6d poly   poly mean %5.2f max %3d   fac mean %6.2f max %4d\n",
            name,
            TIMING ? @sprintf("%8.3f s", a.t) : "",
            a.s.cells,
            sum(a.s.parts),
            a.s.mean_p,
            a.s.max_p,
            a.s.mean_f,
            a.s.max_f
        )
    end
    if TIMING
        @printf(
            "  → approach 2 has %.2fx the cells of approach 1, and is %.2fx %s\n",
            a2.s.cells / a1.s.cells,
            max(a1.t, a2.t) / min(a1.t, a2.t),
            a2.t < a1.t ? "FASTER" : "slower"
        )
    else
        @printf(
            "  → approach 2 has %.2fx the cells of approach 1 (TIMING=0, no times)\n",
            a2.s.cells / a1.s.cells
        )
    end
    return (; label, a1, a2, nodes, tag)
end

primal = run_graph(false)
dual = run_graph(true)

println(
    "\n",
    "="^84,
    "\nRESULT — the cell count reverses, the time verdict does not\n",
    "="^84,
)
@printf(
    "%-22s %10s %10s %11s %10s %10s %9s\n",
    "graph",
    "cells (1)",
    "cells (2)",
    "(2) / (1)",
    "time (1)",
    "time (2)",
    "speed-up"
)
println("-"^84)
for r in (primal, dual)
    tm(x) = TIMING ? @sprintf("%.3fs", x) : "—"
    @printf(
        "%-22s %10d %10d %11.2f %10s %10s %9s\n",
        r.label,
        r.a1.s.cells,
        r.a2.s.cells,
        r.a2.s.cells / r.a1.s.cells,
        tm(r.a1.t),
        tm(r.a2.t),
        TIMING ? @sprintf("%.2fx", r.a1.t / r.a2.t) : "—"
    )
end
println("="^84)
println(
    """
On the PRIMAL graph the node constrains nothing about the future, so every lifted node must still
serve every mode and approach 2 pays |S| copies of the single node's work. It builds more cells and
wins anyway, because on a complete graph the induced common is the union of the node pieces, and
approach 1's cells inherit that union's fragmentation.

On the DUAL graph the node commits to the mode played next, so it refines under one mode instead of
all of them, and approach 2 builds a strictly smaller quotient.""",
)

# ── figures ──────────────────────────────────────────────────────────────────────────────────────
# The same set as experiment 3, deliberately: the two are meant to be read side by side.

# Outside the guard below: `paper_figures_exp2.jl` redraws the same quotients with FIGURES=0, and
# needs the region outlines on them.
function outline!(p)
    for (i, (_, R)) in enumerate(REGIONS)
        plot!(
            p,
            R;
            fillalpha = 0.0,
            linecolor = REGION_COLOURS[i],
            linewidth = 1.3,
            label = "",
        )
    end
    return p
end

if FIGURES
    d = joinpath(@__DIR__, "figures")
    mkpath(d)

    # ACM sigconf puts two columns of 3.33 in on the page, and Plots measures `size` in pixels at
    # 100 dpi, so 333 px is one column. Computer Modern throughout: GR draws maths with its own
    # CM-like engine whatever the text font, so asking for CM everywhere keeps the tick numbers and
    # the labels one family. `paper_figures_exp2.jl` replaces these with acmart's Libertine.
    COL_PX = 333
    paper!(p) = plot!(
        p;
        fontfamily = "Computer Modern",
        tickfontsize = 9,
        guidefontsize = 11,
        legendfontsize = 7,
        grid = false,
        framestyle = :axes,
        left_margin = 1mm,
        right_margin = 1mm,
        top_margin = 1mm,
        bottom_margin = 1mm,
        dpi = 300,
    )

    "Trim the uniform border off a PNG. `pad` leaves a hair of margin so nothing is shaved."
    function crop_png(path; pad = 3)
        img = PNGFiles.load(path)
        bg = img[1, 1]
        rows = findall(i -> any(!=(bg), view(img, i, :)), axes(img, 1))
        cols = findall(j -> any(!=(bg), view(img, :, j)), axes(img, 2))
        (isempty(rows) || isempty(cols)) && return nothing
        r = max(first(rows) - pad, 1):min(last(rows) + pad, size(img, 1))
        c = max(first(cols) - pad, 1):min(last(cols) + pad, size(img, 2))
        PNGFiles.save(path, img[r, c])
        return nothing
    end
    out(fig, nm) =
        (savefig(fig, joinpath(d, nm)); crop_png(joinpath(d, nm)); println("wrote ", nm))

    # Under equal aspect Plots derives the height from the data; fixing it independently aborts.
    function panelsize(p; w = COL_PX, chrome = 58)
        (xlo, xhi), (ylo, yhi) = Plots.xlims(p), Plots.ylims(p)
        return (
            w,
            round(
                Int,
                (w - chrome) * clamp((yhi - ylo) / max(xhi - xlo, eps()), 0.35, 2.2),
            ) + chrome,
        )
    end

    # Approach 1's quotient is a single plane, since the induced common has one node. Set it against
    # the layered figure below.
    for r in (primal, dual)
        p = plot(;
            aspect_ratio = :equal,
            legend = false,
            title = "",
            xlabel = L"x_1",
            ylabel = L"x_2",
        )
        plot!(
            p,
            r.a1.b.quotient;
            what = :states,
            by = :state,
            show_contours = true,
            linewidth = 0.2,
            fillalpha = 0.9,
        )
        outline!(p)
        paper!(p)
        plot!(p; size = panelsize(p))
        out(p, "exp2_$(r.tag)_approach1_quotient.png")
    end

    # Facets per cell, on a log count: most cells are simple, and the construction pays for the tail.
    # `fillrange` pins both series to one baseline, which a log axis otherwise picks per series.
    for r in (primal, dual)
        bins =
            range(0, 1.02 * max(maximum(r.a1.s.faces), maximum(r.a2.s.faces)); length = 45)
        h = histogram(
            r.a1.s.faces;
            bins,
            yscale = :log10,
            fillrange = 0.7,
            alpha = 0.55,
            color = :firebrick,
            linecolor = :firebrick,
            label = "iCLF",
            xlabel = "facets per cell",
            ylabel = "cells",
            legend = :topright,
            title = "",
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
            label = "PCLF",
        )
        vline!(h, [r.a1.s.max_f]; color = :firebrick, ls = :dash, lw = 1.8, label = "")
        vline!(h, [r.a2.s.max_f]; color = :steelblue, ls = :dash, lw = 1.8, label = "")
        paper!(h)
        # wide and short: a histogram carries its information horizontally, and this is the figure
        # that can afford to be short in a two-column layout
        plot!(h; size = (COL_PX, round(Int, 0.44 * COL_PX)))
        out(h, "exp2_$(r.tag)_facets_histogram.png")
    end

    # ── 3-D: approach 2's quotient, the graph node as a vertical axis ────────────────────────────
    # `import`, not `using`: both packages export `plot`, so every Plots figure had to come first,
    # and `using CairoMakie` would make every bare `plot` ambiguous in a REPL holding this file.
    import CairoMakie

    # Makie draws its own text, so these two cannot take the LaTeX route the flat figures use. The
    # next best thing is to hand Makie the same typeface acmart sets: TeX distributions ship Linux
    # Libertine as OTF, which Makie can load, and `kpsewhich` finds it wherever it lives.
    libertine = let kp = Sys.which("kpsewhich")
        p = kp === nothing ? "" : try
                String(strip(read(`$kp LinLibertine_R.otf`, String)))
            catch
                ""
            end
        isfile(p) ? p : nothing
    end

    for r in (primal, dual)
        q = r.a2.b.quotient
        node_z = Dict(nd => z for (nd, z) in zip(r.nodes, (0.0, 1.0)))
        mk = CairoMakie.Figure(;
            size = (900, 700),
            fonts = isnothing(libertine) ? NamedTuple() : (; regular = libertine),
        )
        ax = CairoMakie.Axis3(
            mk[1, 1];
            xlabel = "x₁",
            ylabel = "x₂",
            zlabel = "",
            # the De Bruijn node names say nothing to a reader here, and the axis is only ever two
            # layers, so they are named a and b
            zticks = ([0.0, 1.0], ["node a", "node b"]),
            azimuth = 1.2π,
            elevation = 0.16π,
            xlabelsize = 34,
            ylabelsize = 34,
            zlabelsize = 34,
            xticklabelsize = 28,
            yticklabelsize = 28,
            zticklabelsize = 28,
        )
        DI.plot_augmented_bisimulation!(
            ax,
            q;
            node_z = node_z,
            color_by = :state,
            alpha = 0.35,
            show_contours = false,
        )
        nm = "exp2_$(r.tag)_approach2_3d.png"
        CairoMakie.save(joinpath(d, nm), mk; px_per_unit = 3)
        crop_png(joinpath(d, nm))
        println("wrote ", nm)
    end
end
