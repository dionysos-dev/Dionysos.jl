# EXPERIMENT 3 — what happens when the certificate's node pieces genuinely differ.
#
#     julia --project=. experiment3_diverse_pieces.jl        [RUNGS=n] [SAMPLES=n] [FIGURES=0]
#
# Self-contained: this file defines its own system, certificate, ladder and figures.
#
# Experiments 1 and 2 use certificates whose node pieces are close to one another. This one does
# not, and that is its whole purpose. The working set is the one of experiment 2; only the dynamics
# change, and with them the pieces the solver returns.
#
# The diversity is DISCOVERED, not imposed. Both nodes are given the same conic template of order 2,
# and the solver comes back with two markedly different pieces — a support-function gap of about
# 0.46, against 0.03 for experiment 2's certificate. Forcing the pieces apart by handing each node a
# different template would be the easier route and a much weaker claim.
#
# NO OBSERVATION REGIONS HERE, unlike the other two experiments. With regions, part of the
# refinement serves to respect them, and that part is work both approaches do identically. Removing
# them leaves the certificate's own geometry as the only thing driving the partition, which is
# exactly what this experiment is about. It also removes the specification, so there is nothing to
# solve and nothing to certify: the comparison is on the abstraction alone.
#
# The price is that the ladder has to be given explicitly. The construction stops at the first level
# whose terminal set clears the regions, and with no regions that holds immediately: the ladder would
# be one rung and the quotient one cell per node. So the rung count is fixed, and the outer level is
# taken from a probe that lets approach 1 derive its own, then imposed on approach 2. There is no
# terminal set to overlap anything, so fixing the count is sound here in a way it would not be with
# regions present.
#
# BOTH ORIENTATIONS ARE RUN, and the pair with experiment 2 is the point. Experiment 2 flips the
# orientation with pieces that nearly coincide; this one flips it with pieces that do not. The two
# mechanisms are different and the question is whether approach 2 wins under both:
#
#   primal, complete      the induced common is `min_i V_i`, the UNION of the node pieces, which
#                         fragments as they separate. Approach 1's cells inherit the fragmentation
#                         and approach 2's do not, so approach 2 should win on CELL SIMPLICITY while
#                         building more cells.
#   dual, co-complete     the induced common is `max_i V_i`, the INTERSECTION, which stays convex
#                         however different the pieces are. Approach 1 then pays nothing for the
#                         diversity, so the simplicity channel is closed; what is left is the
#                         per-node reduction, and diversity works against it.

using StaticArrays, LinearAlgebra, JuMP, Clarabel, LazySets, Plots, Printf
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
# The ladder has to be given explicitly here, since the construction's stopping rule is satisfied
# immediately when there is no region to clear. The count is not neutral, and it is chosen on the
# polytope totals rather than on the times, which vary by a factor of three between runs on this
# machine while every count reproduces exactly. Approach 2's saving, as a share of approach 1's
# polytopes, over the two orientations:
#
#   rungs      7        9       10       11
#   primal   46 %     32 %     28 %     23 %
#   dual     -2 %     19 %     21 %     37 %
#
# The two orientations move in opposite directions as the quotient is refined: the primal loses
# margin as the induced common's fragmentation weighs relatively less, the dual gains it as the
# arms' cell counts separate. Seven is where the dual has nothing to show. Nine is the first count
# where both are solidly positive, and it keeps the quotients small enough that the partition
# figures still show individual cells rather than a texture.
const RUNGS = parse(Int, get(ENV, "RUNGS", "9"))
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
# The same working set as experiment 2, deliberately. The two experiments are meant to differ in one
# respect — whether the certificate's node pieces coincide — and giving this one its own geometry
# would add a second difference.
X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
problem = PR.BisimulationQuotientProblem(f, X, typeof(X)[])

println(
    "="^84,
    "\nEXPERIMENT 3 — diverse node pieces, both orientations, no regions\n",
    "="^84,
)

"How different two pieces are: the largest relative gap between their support functions. Exact for
convex sets, needs no volume, and reads as a percentage."
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
    MOI.set(opt, MOI.RawOptimizerAttribute("nb_levels"), RUNGS)
    MOI.set(opt, MOI.RawOptimizerAttribute("max_slices"), RUNGS)
    isnothing(ΓX) || MOI.set(opt, MOI.RawOptimizerAttribute("ΓX"), ΓX)
    MOI.optimize!(opt)
    return (;
        opt,
        quotient = MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
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
    st = PCQ.bisimulation_stats(q)
    return (;
        cells = length(parts),
        parts,
        faces,
        transitions = st[:num_transitions],
        max_p = maximum(parts),
        max_f = maximum(faces),
        mean_p = Statistics.mean(parts),
        mean_f = Statistics.mean(faces),
        sum_f = sum(faces),
    )
end

function run_orientation(dual)
    label = dual ? "DUAL (co-complete)" : "PRIMAL (complete)"
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

    # Approach 1's own natural outer level, then imposed on approach 2. Taking it from approach 2
    # would flatter approach 2; leaving each to its own would have them tile different sets, and the
    # cell counts would measure coverage rather than efficiency.
    ΓX = maximum(MOI.get(build(common).opt, MOI.RawOptimizerAttribute("Γ")))

    println("\n", "#"^84, "\n$label De Bruijn, order 1\n", "#"^84)
    @printf(
        "complete: %-5s  co-complete: %-5s  certified rate %.6f\n",
        PCLF.is_complete(graph, 1:2),
        PCLF.is_co_complete(graph, 1:2),
        pclf.JSRapprox
    )
    @printf(
        "piece gap %.4f   induced common %d convex part(s)   ΓX = %.4f, %d rungs\n",
        gap,
        parts_common,
        ΓX,
        RUNGS
    )
    # The certificate itself. A polyhedral piece is the gauge `V_s(x) = max_i |(G x)_i| / w_i`, so
    # its Γ-sublevel set is the symmetric polytope `{x : |G x| ≤ Γ w}`.
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

    # Build both arms first, then time both. A timed region pays garbage collection proportional to
    # the live heap, so an arm timed while the other does not yet exist is timed on a lighter heap:
    # building and timing in one pass favours whichever arm goes first, by enough to reverse a
    # verdict. Both quotients are alive for both measurements here.
    arm(cert) = (; b = build(cert; ΓX = ΓX), cert)
    a1, a2 = arm(common), arm(pclf)
    a1 = (;
        a1...,
        s = stats(a1.b.quotient),
        t = timed(() -> build(a1.cert; ΓX = ΓX), SAMPLES),
    )
    a2 = (;
        a2...,
        s = stats(a2.b.quotient),
        t = timed(() -> build(a2.cert; ΓX = ΓX), SAMPLES),
    )
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
        a2.t < a1.t ? "FASTER" : "SLOWER"
    )
    return (; label, tag = dual ? "dual" : "primal", nodes, gap, parts_common, a1, a2)
end

primal = run_orientation(false)
dual = run_orientation(true)

println("\n", "="^84, "\nRESULT\n", "="^84)
@printf(
    "%-20s %6s %6s %9s %9s %8s %8s %8s %9s\n",
    "orientation",
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
        "%-20s %6.3f %6d %9d %9d %8.2f %8d %8d %8.2fx\n",
        r.label,
        r.gap,
        r.parts_common,
        r.a1.s.cells,
        r.a2.s.cells,
        r.a2.s.cells / r.a1.s.cells,
        r.a1.s.max_p,
        r.a2.s.max_p,
        r.a1.t / r.a2.t
    )
end
println("="^84)
println(
    """
Read this beside experiment 2, which flips the same orientation with pieces that nearly coincide.
The two channels behind approach 2's advantage are separate, and the certificate decides which of
them is open.

On the PRIMAL graph the induced common is the union of the node pieces, so it fragments as they
separate: approach 1's cells inherit the fragmentation and approach 2's do not. That channel is wide
open here, and it is what the `maxp` columns measure.

On the DUAL graph the induced common is their intersection, which stays convex however different the
pieces are. Approach 1 pays nothing for the diversity, the simplicity channel is closed, and all that
is left is the per-node reduction — which diversity works against. Whether approach 2 still comes out
ahead there is what the speed-up column answers.""",
)

if FIGURES
    d = joinpath(@__DIR__, "figures")
    mkpath(d)

    # Under equal aspect Plots derives the height from the width and the data, so a height fixed
    # independently of the data over-constrains the layout and Plots aborts.
    function panelsize(p; w = 520, chrome = 95)
        (xlo, xhi), (ylo, yhi) = Plots.xlims(p), Plots.ylims(p)
        return (
            w,
            round(Int, w * clamp((yhi - ylo) / max(xhi - xlo, eps()), 0.35, 2.2)) + chrome,
        )
    end

    # Approach 1's quotient: one plane, because the induced common has a single node. This is the
    # object the earlier method works on, and the figure to set against the layered one beside it.
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
        plot!(p; size = panelsize(p))
        savefig(p, joinpath(d, "exp3_$(r.tag)_approach1_quotient.png"))
        println("wrote exp3_$(r.tag)_approach1_quotient.png")
    end

    # Facets per cell. Log counts: most cells are simple, the tail is not, and the tail is what the
    # construction pays for. Drawn before the 3-D block, since `yscale = :log10` has been observed to
    # disturb later Plots figures and everything Plots-based must precede the CairoMakie import anyway.
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
        savefig(h, joinpath(d, "exp3_$(r.tag)_facets_histogram.png"))
        println("wrote exp3_$(r.tag)_facets_histogram.png")
    end

    # ── 3-D: approach 2's quotient, with the graph node as a vertical axis ───────────────────────
    # Both packages export `plot`, so every Plots figure above had to come first — and the import
    # is deliberate: `using CairoMakie` makes every bare `plot` in this session ambiguous, which
    # breaks redrawing a figure from a REPL that has the script loaded.
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
        CairoMakie.save(joinpath(d, "exp3_$(r.tag)_approach2_3d.png"), mk; px_per_unit = 3)
        println("wrote exp3_$(r.tag)_approach2_3d.png")
    end
end
