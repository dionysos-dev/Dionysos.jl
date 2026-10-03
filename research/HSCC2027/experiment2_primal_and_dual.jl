# EXPERIMENT 2 — the same system on the primal and on the dual De Bruijn graph.
#
#     julia --project=. experiment2_primal_and_dual.jl        [SAMPLES=n] [FIGURES=0]
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

using StaticArrays, JuMP, Clarabel, LazySets, Plots, Printf
import HybridSystems, Statistics
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
        b = build(cert, NB, ΓX)
        parts = b.D isa LazySets.UnionSetArray ? b.D.array : [UT._as_hpolytope(b.D)]
        for (nm, R) in REGIONS, P in parts
            LazySets.isdisjoint(P, UT._as_hpolytope(R)) ||
                error("[$label / $name] the terminal set D meets region $nm")
        end
        return (; b, s = stats(b.quotient), cert)
    end
    a1, a2 = arm("approach 1", common), arm("approach 2", pclf)
    t1 = timed(() -> build(a1.cert, NB, ΓX), SAMPLES)
    t2 = timed(() -> build(a2.cert, NB, ΓX), SAMPLES)
    a1, a2 = (; a1..., t = t1), (; a2..., t = t2)
    for (name, a) in (("approach 1", a1), ("approach 2", a2))
        @printf(
            "  %-10s %8.3f s  %5d cells  %6d poly   poly mean %5.2f max %3d   fac mean %6.2f max %4d\n",
            name,
            a.t,
            a.s.cells,
            sum(a.s.parts),
            a.s.mean_p,
            a.s.max_p,
            a.s.mean_f,
            a.s.max_f
        )
    end
    @printf(
        "  → approach 2 has %.2fx the cells of approach 1, and is %.2fx %s\n",
        a2.s.cells / a1.s.cells,
        max(a1.t, a2.t) / min(a1.t, a2.t),
        a2.t < a1.t ? "FASTER" : "slower"
    )
    return (; label, a1, a2, nodes, tag = dual ? "dual" : "primal")
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
    @printf(
        "%-22s %10d %10d %11.2f %9.3fs %9.3fs %8.2fx\n",
        r.label,
        r.a1.s.cells,
        r.a2.s.cells,
        r.a2.s.cells / r.a1.s.cells,
        r.a1.t,
        r.a2.t,
        r.a1.t / r.a2.t
    )
end
println("="^84)
println(
    """
On the PRIMAL graph the node constrains nothing about the future, so every lifted node must still
serve every mode and approach 2 pays |S| copies of the single node's work. It builds more cells and
wins anyway, because determinising a complete graph produces a non-convex union and approach 1's
cells inherit its fragmentation.

On the DUAL graph the node commits to the mode played next, so it refines under one mode instead of
all of them, and approach 2 builds a strictly smaller quotient.""",
)

# ── figures ──────────────────────────────────────────────────────────────────────────────────────
# The same set as experiment 3, deliberately: the two are meant to be read side by side.
if FIGURES
    d = joinpath(@__DIR__, "figures")
    mkpath(d)

    function outline!(p)
        for (i, (_, R)) in enumerate(REGIONS)
            plot!(
                p,
                R;
                fillalpha = 0.0,
                linecolor = REGION_COLOURS[i],
                linewidth = 2.0,
                label = "",
            )
        end
        return p
    end
    # Under equal aspect Plots derives the height from the data; fixing it independently aborts.
    function panelsize(p; w = 520, chrome = 95)
        (xlo, xhi), (ylo, yhi) = Plots.xlims(p), Plots.ylims(p)
        return (
            w,
            round(Int, w * clamp((yhi - ylo) / max(xhi - xlo, eps()), 0.35, 2.2)) + chrome,
        )
    end

    # Approach 1's quotient is a single plane, since the induced common has one node. Set it against
    # the layered figure below.
    for r in (primal, dual)
        p = plot(;
            aspect_ratio = :equal,
            legend = false,
            titlefontsize = 10,
            title = "$(r.label) — induced common, $(r.a1.s.cells) cells, " *
                    "worst $(r.a1.s.max_p) parts",
        )
        plot!(
            p,
            r.a1.b.quotient;
            what = :states,
            by = :state,
            show_contours = true,
            linewidth = 0.3,
            fillalpha = 0.9,
        )
        outline!(p)
        plot!(p; size = panelsize(p))
        savefig(p, joinpath(d, "exp2_$(r.tag)_approach1_quotient.png"))
        println("wrote exp2_$(r.tag)_approach1_quotient.png")
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
            label = "approach 1 (determinise first)",
            xlabel = "facets per cell",
            ylabel = "cells",
            legend = :topright,
            framestyle = :box,
            grid = false,
            ylims = (0.7, 5e3),
            size = (660, 410),
            title = r.label,
            titlefontsize = 10,
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
            label = "approach 2 (on the PCLF)",
        )
        vline!(h, [r.a1.s.max_f]; color = :firebrick, ls = :dash, lw = 1.5, label = "")
        vline!(h, [r.a2.s.max_f]; color = :steelblue, ls = :dash, lw = 1.5, label = "")
        savefig(h, joinpath(d, "exp2_$(r.tag)_facets_histogram.png"))
        println("wrote exp2_$(r.tag)_facets_histogram.png")
    end

    # ── 3-D: approach 2's quotient, the graph node as a vertical axis ────────────────────────────
    # `import`, not `using`: both packages export `plot`, so every Plots figure had to come first,
    # and `using CairoMakie` would make every bare `plot` ambiguous in a REPL holding this file.
    import CairoMakie

    for r in (primal, dual)
        q = r.a2.b.quotient
        node_z = Dict(nd => z for (nd, z) in zip(r.nodes, (0.0, 1.0)))
        mk = CairoMakie.Figure(; size = (900, 700))
        ax = CairoMakie.Axis3(
            mk[1, 1];
            xlabel = "x₁",
            ylabel = "x₂",
            zlabel = "memory (graph node)",
            zticks = ([0.0, 1.0], ["node $(r.nodes[1])", "node $(r.nodes[2])"]),
            azimuth = 1.2π,
            elevation = 0.16π,
            title = "$(r.label) — PCLF quotient, $(r.a2.s.cells) cells, " *
                    "worst $(r.a2.s.max_p) parts",
        )
        DI.plot_augmented_bisimulation!(
            ax,
            q;
            node_z = node_z,
            color_by = :state,
            alpha = 0.35,
            show_contours = false,
        )
        CairoMakie.save(joinpath(d, "exp2_$(r.tag)_approach2_3d.png"), mk; px_per_unit = 3)
        println("wrote exp2_$(r.tag)_approach2_3d.png")
    end
end
