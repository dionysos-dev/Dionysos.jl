# EXPERIMENT 2 — the same system on both orientations of the De Bruijn graph.
#
#     julia --project=. experiment2_graph_orientation.jl        [SAMPLES=n] [FIGURES=0]
#
# Self-contained: this file defines its own system, certificate and ladder. It shares nothing with
# experiment 1, which needs different choices — in particular a different rule for the outer level,
# because there the observation regions are larger than the nominal working set and here they are
# much smaller.
#
# One system, one certificate family, one working set. The only thing that changes is what a graph
# node means:
#
#   primal   the node records the mode just PLAYED      — every node emits every mode (complete)
#   dual     the node commits to the mode played NEXT   — every node emits one mode (co-complete)
#
# Approach 2 lifts the abstraction over the graph's nodes, so the obvious expectation is that it
# always builds a bigger quotient. It does not. Flipping the orientation reverses the cell count
# while leaving the time verdict alone:
#
#   primal   approach 2 builds MORE cells than approach 1 and is still faster
#   dual     approach 2 builds FEWER cells than approach 1 and is faster
#
# That is the point of the experiment: the number of cells does not predict the cost. What governs
# the reversal, and how to read it off the graph, is left to future work.
#
# A note on the certificate, because the choice is not free. On a complete graph the induced common
# is `min_i V_i`, whose sublevel set is the UNION of the node pieces — non-convex, and it fragments
# as the pieces separate. On a co-complete graph it is `max_i V_i`, the INTERSECTION, which stays
# convex however different the pieces are. Piece diversity therefore helps approach 2 on the primal
# graph and hurts it on the dual: forcing the pieces apart on a dual graph has been measured to put
# approach 2's cell ratio at 0.29–0.41, that is, approach 2 building two to three times MORE cells.
# The template here is shared between the nodes, which keeps the pieces close and lets the dual show
# the effect it is here to show.

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
const ATOL = 1e-3
const THETA = π / 6        # the shared template's orientation
const SOLVER =
    JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000, "verbose" => false)
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

println("="^84, "\nEXPERIMENT 2 — one system, both graph orientations\n", "="^84)

# ── the two approaches, on one orientation ───────────────────────────────────────────────────────
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
        cells = length(parts), parts, faces,
        max_p = maximum(parts), max_f = maximum(faces),
        mean_p = Statistics.mean(parts), mean_f = Statistics.mean(faces),
    )
end

function run_orientation(dual)
    label = dual ? "DUAL (co-complete)" : "PRIMAL (complete)"
    graph = PCLF.generate_DeBruijn_edges(2, 1; dual = dual)
    nodes = sort(collect(graph.verts); by = string)
    Rot = [cos(THETA) -sin(THETA); sin(THETA) cos(THETA)]
    pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
        f, graph, SOLVER; Gmats = Dict(nd => Rot for nd in nodes), MLF = true, verbose = false,
    )
    isinf(pclf.JSRapprox) && error("no certificate on the $label graph")
    common = PCLF.build_common_lyapunov(pclf)
    γ = pclf.JSRapprox

    # The level at which the induced COMMON covers the working set, imposed on both approaches. The
    # certificate is left at its natural scale, whose outer level is already of order 10^3: the level
    # tolerances handle that, and rescaling would change the row norms and hence the geometric
    # thickness of `atol`.
    ΓX = PCQ.compute_tau_X(common, UT._as_hpolytope(X))
    # The construction's own stopping rule: the first level whose terminal set no longer meets a
    # region, for both certificates.
    NB = something(
        findfirst(
            j ->
                PCQ.all_nodes_clear_regions(pclf, ΓX * γ^(j - 1), REGIONS_H; tol = 1e-2) &&
                PCQ.all_nodes_clear_regions(common, ΓX * γ^(j - 1), REGIONS_H; tol = 1e-2),
            1:60,
        ),
        -1,
    )
    NB > 0 || error("no region-free terminal level on the $label graph")

    common_parts = let S = PCLF.get_sublevel_set(common.pieces[:clf], 1.0; atol = 1e-6)
        S isa LazySets.UnionSetArray ? length(S.array) : 1
    end

    println("\n", "#"^84, "\n$label De Bruijn, order 1\n", "#"^84)
    @printf("complete: %-5s  co-complete: %-5s  certified rate %.6f\n",
            PCLF.is_complete(graph, 1:2), PCLF.is_co_complete(graph, 1:2), γ)
    @printf("induced common: %d convex part(s)   ladder: ΓX = %.4f, %d slices\n",
            common_parts, ΓX, NB)

    # Build both arms first, then time both. A timed region pays garbage collection proportional to
    # the live heap, so an arm timed while the other does not yet exist is timed on a lighter heap:
    # building and timing in one pass favours whichever arm goes first, by enough to reverse a
    # verdict. Both quotients are alive for both measurements here.
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
        @printf("  %-10s %8.3f s  %5d cells  %6d poly   poly mean %5.2f max %3d   fac mean %6.2f max %4d\n",
                name, a.t, a.s.cells, sum(a.s.parts), a.s.mean_p, a.s.max_p, a.s.mean_f, a.s.max_f)
    end
    @printf("  → approach 2 has %.2fx the cells of approach 1, and is %.2fx %s\n",
            a2.s.cells / a1.s.cells, max(a1.t, a2.t) / min(a1.t, a2.t),
            a2.t < a1.t ? "FASTER" : "slower")
    return (; label, a1, a2, nodes, tag = dual ? "dual" : "primal")
end

primal = run_orientation(false)
dual = run_orientation(true)

println("\n", "="^84, "\nRESULT — the cell count reverses, the time verdict does not\n", "="^84)
@printf("%-22s %10s %10s %11s %10s %10s %9s\n",
        "orientation", "cells (1)", "cells (2)", "(2) / (1)", "time (1)", "time (2)", "speed-up")
println("-"^84)
for r in (primal, dual)
    @printf("%-22s %10d %10d %11.2f %9.3fs %9.3fs %8.2fx\n",
            r.label, r.a1.s.cells, r.a2.s.cells, r.a2.s.cells / r.a1.s.cells,
            r.a1.t, r.a2.t, r.a1.t / r.a2.t)
end
println("="^84)
println("""
On the PRIMAL graph the node constrains nothing about the future, so every lifted node must still
serve every mode and approach 2 pays |S| copies of the single node's work. It builds more cells and
wins anyway, because determinising a complete graph produces a non-convex union and approach 1's
cells inherit its fragmentation.

On the DUAL graph the node commits to the mode played next, so it refines under one mode instead of
all of them, and approach 2 builds a strictly smaller quotient.""")

# ── figures ──────────────────────────────────────────────────────────────────────────────────────
# Deliberately the same set as experiment 3, which differs from this one only in how far apart the
# node pieces sit: the two are meant to be read side by side.
if FIGURES
    d = joinpath(@__DIR__, "figures")
    mkpath(d)

    function outline!(p)
        for (i, (_, R)) in enumerate(REGIONS)
            plot!(p, R; fillalpha = 0.0, linecolor = REGION_COLOURS[i], linewidth = 2.0, label = "")
        end
        return p
    end
    # Under equal aspect Plots derives the height from the width and the data, so a height fixed
    # independently of the data over-constrains the layout and Plots aborts.
    function panelsize(p; w = 520, chrome = 95)
        (xlo, xhi), (ylo, yhi) = Plots.xlims(p), Plots.ylims(p)
        return (w, round(Int, w * clamp((yhi - ylo) / max(xhi - xlo, eps()), 0.35, 2.2)) + chrome)
    end

    # Approach 1's quotient: one plane, because the induced common has a single node. This is the
    # object the earlier method works on, and the figure to set against the layered one below it.
    for r in (primal, dual)
        p = plot(; aspect_ratio = :equal, legend = false, titlefontsize = 10,
                 title = "$(r.label) — induced common, $(r.a1.s.cells) cells, " *
                         "worst $(r.a1.s.max_p) parts")
        plot!(p, r.a1.b.quotient; what = :states, by = :state, show_contours = true,
              linewidth = 0.3, fillalpha = 0.9, merge_series = false)
        outline!(p)
        plot!(p; size = panelsize(p))
        savefig(p, joinpath(d, "exp2_$(r.tag)_approach1_quotient.png"))
        println("wrote exp2_$(r.tag)_approach1_quotient.png")
    end

    # Facets per cell. Log counts: most cells are simple, the tail is not, and the tail is what the
    # construction pays for. Drawn before the 3-D block, since `yscale = :log10` has been observed to
    # disturb later Plots figures and everything Plots-based must precede the CairoMakie import anyway.
    for r in (primal, dual)
        bins = range(0, 1.02 * max(maximum(r.a1.s.faces), maximum(r.a2.s.faces)); length = 45)
        h = histogram(r.a1.s.faces; bins, yscale = :log10, fillrange = 0.7, alpha = 0.55,
                      color = :firebrick, linecolor = :firebrick,
                      label = "approach 1 (determinise first)",
                      xlabel = "facets per cell", ylabel = "cells", legend = :topright,
                      framestyle = :box, grid = false, ylims = (0.7, 5e3), size = (660, 410),
                      title = r.label, titlefontsize = 10)
        histogram!(h, r.a2.s.faces; bins, yscale = :log10, fillrange = 0.7, alpha = 0.55,
                   color = :steelblue, linecolor = :steelblue, label = "approach 2 (on the PCLF)")
        vline!(h, [r.a1.s.max_f]; color = :firebrick, ls = :dash, lw = 1.5, label = "")
        vline!(h, [r.a2.s.max_f]; color = :steelblue, ls = :dash, lw = 1.5, label = "")
        savefig(h, joinpath(d, "exp2_$(r.tag)_facets_histogram.png"))
        println("wrote exp2_$(r.tag)_facets_histogram.png")
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
            mk[1, 1]; xlabel = "x₁", ylabel = "x₂", zlabel = "memory (graph node)",
            zticks = ([0.0, 1.0], ["node $(r.nodes[1])", "node $(r.nodes[2])"]),
            azimuth = 1.2π, elevation = 0.16π,
            title = "$(r.label) — PCLF quotient, $(r.a2.s.cells) cells, " *
                    "worst $(r.a2.s.max_p) parts",
        )
        DI.plot_augmented_bisimulation!(ax, q; node_z = node_z, color_by = :state,
                                        alpha = 0.35, show_contours = false)
        CairoMakie.save(joinpath(d, "exp2_$(r.tag)_approach2_3d.png"), mk; px_per_unit = 3)
        println("wrote exp2_$(r.tag)_approach2_3d.png")
    end
end
