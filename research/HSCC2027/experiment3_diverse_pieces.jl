# EXPERIMENT 3 — the same two graphs, with node pieces that do not coincide.
#
#     julia --project=. experiment3_diverse_pieces.jl        [LEVELS=n] [SAMPLES=n] [FIGURES=0]
#
# Self-contained: this file defines its own system, certificate, levels and figures.
#
# Experiments 1 and 2 use certificates whose node pieces are close to one another; this one changes
# the dynamics so that they are not. The working set is experiment 2's, deliberately, so that the
# pieces are the only difference between the two. The diversity is DISCOVERED, not imposed: both
# nodes get the same conic template of order 2 and the solver returns two markedly different pieces.
#
# NO OBSERVATION REGIONS HERE. With regions, part of the refinement serves to respect them, and that
# is work both approaches do identically; removing them leaves the certificate's own geometry as the
# only thing driving the partition. It also removes the specification, so the comparison is on the
# abstraction alone, and the stopping rule, so the depth is given by LEVELS instead of derived.
#
# Both graphs are run, and reading this beside experiment 2 separates the two channels behind
# approach 2's advantage: on the primal graph the induced common is the UNION of the node pieces,
# which fragments as they separate, so approach 2 wins on cell simplicity; on the dual it is their
# INTERSECTION, which stays convex, so that channel is closed and only the per-node reduction is
# left.

using JuMP, Clarabel, LazySets, Plots, Printf, LaTeXStrings
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
# Each graph builds seven quotients: one to fix the outer level, one per arm, and a warm-up plus
# SAMPLES timed rounds per arm. CACHE=1 keeps the first three on disk, which is what makes redrawing
# the figures cheap; TIMING=0 drops the other four, which are the measurement and are never cached.
const USE_CACHE = get(ENV, "CACHE", "1") == "1"
const TIMING = get(ENV, "TIMING", "1") == "1"
const LEVELS = parse(Int, get(ENV, "LEVELS", "9"))
const ATOL = 1e-3
const SOLVER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)

# ── the problem ──────────────────────────────────────────────────────────────────────────────────
A1 = (1.0 / 10.0) * [1.5519 0.4474; 7.6412 7.4716]
A2 = (1.0 / 10.0) * [0.4750 9.1755; 1.8955 0.1850]
f = ST.with_switching(
    HybridSystems.discreteswitchedsystem([A1, A2]),
    HybridSystems.ControlledSwitching(),
)
X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
problem = PR.BisimulationQuotientProblem(f, X, typeof(X)[])

println(
    "="^84,
    "\nEXPERIMENT 3 — diverse node pieces, the primal and the dual graph, no regions\n",
    "="^84,
)

"How different two pieces are: the largest relative gap between their support functions. Exact for
convex sets and needs no volume."
function piece_gap(P1, P2; K = 64)
    return maximum(
        abs(LazySets.ρ(d, P1) - LazySets.ρ(d, P2)) /
        max(abs(LazySets.ρ(d, P1)), abs(LazySets.ρ(d, P2)), eps()) for
        d in ([cos(2π * k / K), sin(2π * k / K)] for k in 0:(K - 1))
    )
end

function build(cert; ΓX = nothing)
    opt = MOI.instantiate(PCQ.OptimizerBisimulationQuotient)
    MOI.set(opt, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("pclf"), cert)
    MOI.set(opt, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.set(opt, MOI.RawOptimizerAttribute("atol"), ATOL)
    MOI.set(opt, MOI.RawOptimizerAttribute("nb_levels"), LEVELS)
    MOI.set(opt, MOI.RawOptimizerAttribute("max_slices"), LEVELS)
    isnothing(ΓX) || MOI.set(opt, MOI.RawOptimizerAttribute("ΓX"), ΓX)
    MOI.optimize!(opt)
    return (;
        opt,
        quotient = MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
    )
end

"The quotient `cert` builds, from cache/ when it is there. Timed builds never come through here."
function build_cached(name, tag, cert, ΓX)
    path = joinpath(@__DIR__, "cache", "exp3_$(tag)_$(name)_$(LEVELS).jld2")
    if USE_CACHE && isfile(path)
        return JLD2.jldopen(path, "r") do file
            return (; quotient = file["quotient"])
        end
    end
    b = build(cert; ΓX = ΓX)
    if USE_CACHE
        mkpath(dirname(path))
        JLD2.jldopen(path, "w") do file
            return file["quotient"] = b.quotient
        end
    end
    return (; b.quotient)
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
    pclf = PCLF.compute_polyhedral_pieces_pclf(
        f,
        graph,
        SOLVER,
        PCLF.conic_partitions_dict_2d(2, nodes);
        MLF = true,
    )
    isinf(pclf.JSRapprox) && error("no certificate on the $label graph")
    common = PCLF.build_common_lyapunov(pclf)

    parts_common = let S = PCLF.get_sublevel_set(common.pieces[:clf], 1.0; atol = 1e-6)
        S isa LazySets.UnionSetArray ? length(S.array) : 1
    end
    gap = piece_gap(
        PCLF.get_sublevel_set(pclf.pieces[nodes[1]], 1.0),
        PCLF.get_sublevel_set(pclf.pieces[nodes[2]], 1.0),
    )

    # Approach 1's own natural outer level, imposed on approach 2. Taking it from approach 2 would
    # flatter approach 2, and leaving each to its own would have them tile different sets.
    ΓX = let path = joinpath(@__DIR__, "cache", "exp3_$(tag)_GX_$(LEVELS).jld2")
        if USE_CACHE && isfile(path)
            JLD2.jldopen(path, "r") do file
                return file["GX"]
            end
        else
            v = maximum(MOI.get(build(common).opt, MOI.RawOptimizerAttribute("Γ")))
            if USE_CACHE
                mkpath(dirname(path))
                JLD2.jldopen(path, "w") do file
                    return file["GX"] = v
                end
            end
            v
        end
    end

    println("\n", "#"^84, "\n$label De Bruijn, order 1\n", "#"^84)
    @printf(
        "complete: %-5s  co-complete: %-5s  certified rate %.6f\n",
        PCLF.is_complete(graph, 1:2),
        PCLF.is_co_complete(graph, 1:2),
        pclf.JSRapprox
    )
    @printf(
        "piece gap %.4f   induced common %d convex part(s)   ΓX = %.4f, %d levels\n",
        gap,
        parts_common,
        ΓX,
        LEVELS
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
    arm(name, cert) = (; b = build_cached(name, tag, cert, ΓX), cert)
    a1, a2 = arm("approach1", common), arm("approach2", pclf)
    time_arm(a) = TIMING ? timed(() -> build(a.cert; ΓX = ΓX), SAMPLES) : NaN
    a1 = (; a1..., s = stats(a1.b.quotient), t = time_arm(a1))
    a2 = (; a2..., s = stats(a2.b.quotient), t = time_arm(a2))
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
            a2.t < a1.t ? "FASTER" : "SLOWER"
        )
    else
        @printf(
            "  → approach 2 has %.2fx the cells of approach 1 (TIMING=0, no times)\n",
            a2.s.cells / a1.s.cells
        )
    end
    return (; label, tag, nodes, gap, parts_common, a1, a2)
end

primal = run_graph(false)
dual = run_graph(true)

println("\n", "="^84, "\nRESULT\n", "="^84)
@printf(
    "%-20s %6s %6s %9s %9s %8s %8s %8s %9s\n",
    "graph",
    "gap",
    "parts",
    "cells (1)",
    "cells (2)",
    "(2)/(1)",
    "maxp (1)",
    "maxp (2)",
    "speed-up"
)
println("-"^84)
for r in (primal, dual)
    @printf(
        "%-20s %6.3f %6d %9d %9d %8.2f %8d %8d %9s\n",
        r.label,
        r.gap,
        r.parts_common,
        r.a1.s.cells,
        r.a2.s.cells,
        r.a2.s.cells / r.a1.s.cells,
        r.a1.s.max_p,
        r.a2.s.max_p,
        TIMING ? @sprintf("%.2fx", r.a1.t / r.a2.t) : "—"
    )
end
println("="^84)
println(
    """
Read this beside experiment 2, which runs the same two graphs with pieces that nearly coincide. The
two channels behind approach 2's advantage are separate, and the certificate decides which is open.

On the PRIMAL graph the induced common is the union of the node pieces, so it fragments as they
separate: approach 1's cells inherit the fragmentation and approach 2's do not. That is what the
`maxp` columns measure.

On the DUAL graph the induced common is their intersection, which stays convex however different the
pieces are. Approach 1 pays nothing for the diversity, so the simplicity channel is closed and only
the per-node reduction is left — which diversity works against. Whether approach 2 still comes out
ahead is what the speed-up column answers.""",
)

# ── figures ──────────────────────────────────────────────────────────────────────────────────────
# The same set as experiment 2, deliberately: the two are meant to be read side by side.
if FIGURES
    d = joinpath(@__DIR__, "figures")
    mkpath(d)

    # ACM sigconf puts two columns of 3.33 in on the page, and Plots measures `size` in pixels at
    # 100 dpi, so 333 px is one column. Computer Modern throughout: GR draws maths with its own
    # CM-like engine whatever the text font, so asking for CM everywhere keeps the tick numbers and
    # the labels one family. `paper_figures_exp3.jl` replaces these with acmart's Libertine.
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
        paper!(p)
        plot!(p; size = panelsize(p))
        out(p, "exp3_$(r.tag)_approach1_quotient.png")
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
        out(h, "exp3_$(r.tag)_facets_histogram.png")
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
        nm = "exp3_$(r.tag)_approach2_3d.png"
        CairoMakie.save(joinpath(d, nm), mk; px_per_unit = 3)
        crop_png(joinpath(d, nm))
        println("wrote ", nm)
    end
end
