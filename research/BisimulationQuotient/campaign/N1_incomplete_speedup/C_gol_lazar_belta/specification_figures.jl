# GLB primal — the certified region of their co-safe formula, both arms, both quantifiers.
#
#     julia --project=test .../C_gol_lazar_belta/specification_figures.jl
#
# Builds are CACHED next to this script (`cache_route{1,2}.jld2`): route 1 takes about 26 minutes and
# nothing below changes it, so the first run pays and every later one starts from the cache. Delete
# the files to force a rebuild.
#
# TWO QUESTIONS THIS ANSWERS
#
# 1. **Why is the specification solved faster on the LARGER quotient?** Route 2 carries 10 158 cells
#    against route 1's 6 613 and still solves ∃ and ∀ in half the time. The fixed point itself is
#    combinatorial and route 2's is the bigger one — so the answer is not in the fixed point. It is
#    in the operations that touch geometry, and those are charged per *convex part*, not per cell:
#    locating the initial set tests every cell for disjointness, and the closed-loop simulation
#    locates a point in the partition at every step. The script prints Σ parts and Σ facets for both
#    arms, which is the quantity those operations actually pay.
#
# 2. **What does the certified region look like?** Green where the specification is certified, red
#    where it is not. Route 1 has one plane; route 2 has one per graph node, which is what the 3-D
#    views show with the node as a vertical axis — a concrete state sits on two layers at once and
#    the layers need not agree.
#
# `print_level = 1` on the solves is deliberate. On a graph whose nodes enable different modes the
# universal solve completes the missing moves with sink states, and a silent completion once cost
# this project a whole experiment (see `report_completion` in `cosafe_ltl_problem.jl`). The primal
# De Bruijn graph is complete, so the report should show **no completions** — and seeing that is the
# point of asking for it.

include(joinpath(dirname(dirname(dirname(@__DIR__))), "common.jl"))

using Printf
using Spot
import Statistics

gr()

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)
const NB_LEVELS = 7          # the first rung count whose terminal set clears R1, R2, R3
const GAMMA_X = 6.1663       # route 1's own natural cover of their X, imposed on both arms

(; f, R1, R2, R3) = gol_lazar_belta_problem()
X = let L = gol_lazar_belta_L()
    LazySets.HPolytope(vcat(L, -L), fill(10.0, 2 * size(L, 1)))   # their Γ_X = 10
end
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])
φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)     # point `a` of the paper

println("computing the PCLF (the oracle's given) ...")
graph = PCLF.generate_DeBruijn_edges(2, 1; dual = false)
nodes = sort(collect(graph.verts); by = string)
pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    graph,
    OPTIMIZER,
    PCLF.conic_partitions_dict_2d(2, nodes);
    MLF = true,
)
common = PCLF.build_common_lyapunov(pclf)
@printf("  rate %.6f, %d nodes\n\n", pclf.JSRapprox, length(nodes))

function cached(name, cert)
    path = joinpath(@__DIR__, "cache_$name.jld2")
    if isfile(path)
        println("loading $name from cache ...")
        opt = import_optimizer_jld2(path)
        return (;
            quotient = MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
            D = MOI.get(opt, MOI.RawOptimizerAttribute("D")),
        )
    end
    println("building $name (not cached — route 1 takes ~26 min) ...")
    built = build_quotient(
        problem,
        cert;
        atol = 1e-3,
        nb_levels = NB_LEVELS,
        ΓX = GAMMA_X,
        max_slices = NB_LEVELS,
        print_level = 0,
    )
    export_optimizer_jld2(built.optimizer, path)
    return built
end

r1 = cached("route1", common)
r2 = cached("route2", pclf)

# ── QUESTION 1: where the solve's time actually goes ──────────────────────────────────────────────
# Cells are the unit of the fixed point; convex PARTS are the unit of every geometric test. If the
# larger quotient carries fewer parts in total, "more cells" and "faster solve" stop being in
# tension.
println("\n", "="^78)
println("CELL COUNT IS NOT WHAT THE SOLVE PAYS FOR")
println("="^78)
@printf(
    "%-10s %8s %10s %10s %9s %9s\n",
    "arm",
    "cells",
    "Σ parts",
    "Σ facets",
    "mean p/c",
    "max p/c"
)
println("-"^78)
for (nm, b) in (("route 1", r1), ("route 2", r2))
    parts, faces = PCQ.cell_complexities(b.quotient)
    @printf(
        "%-10s %8d %10d %10d %9.2f %9d\n",
        nm,
        length(parts),
        sum(parts),
        sum(faces),
        Statistics.mean(parts),
        maximum(parts)
    )
end
println("="^78)

# ── The specification, both arms, both quantifiers ───────────────────────────────────────────────
ap_to_obs = Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3)
f_forall = ST.with_switching(f, HybridSystems.AutonomousSwitching())

function solve(f_sys, b, label)
    println("\n--- $label ---")
    t0 = time()
    res = synthesize_cosafe_ltl(
        f_sys,
        b.quotient,
        Dionysos.spot_stepper(φ),
        Dict{Symbol, Any}(:D => b.D, :R1 => R1, :R2 => R2, :R3 => R3),
        ap_to_obs,
        x0;
        print_level = 1,
    )
    @printf(
        "  %.3f s, %d of %d cells certified\n",
        time() - t0,
        length(res.controllable_set),
        length(b.quotient.states)
    )
    return res
end

syn1 = solve(f, r1, "route 1 — synthesis (∃)")
ver1 = solve(f_forall, r1, "route 1 — verification (∀)")
syn2 = solve(f, r2, "route 2 — synthesis (∃)")
ver2 = solve(f_forall, r2, "route 2 — verification (∀)")

# ── Planar figures ───────────────────────────────────────────────────────────────────────────────
save(fig, name) =
    (savefig(fig, joinpath(@__DIR__, name * ".png")); println("wrote $name.png"))

# Route 1 lives on one plane: its induced common has a single node.
save(
    plot_synthesis_result(
        r1.quotient,
        syn1,
        problem,
        "Route 1 (induced common) — synthesis ∃",
    ),
    "spec_route1_synthesis",
)
save(
    plot_synthesis_result(
        r1.quotient,
        ver1,
        problem,
        "Route 1 (induced common) — verification ∀",
    ),
    "spec_route1_verification",
)

# Route 2 carries one partition per node. Merged they would superimpose two planes and show neither,
# so each node gets its own panel: a concrete state sits on both, and the two need not agree.
function node_panels(q, res, ttl)
    panels = map(quotient_nodes(q)) do nd
        ids(set) = [i for i in set if q.states[i].node == nd]
        p = plot(;
            aspect_ratio = :equal,
            legend = false,
            title = "$ttl — node $nd",
            titlefontsize = 9,
        )
        _plot_winning_losing!(
            p,
            q,
            ids(res.controllable_set),
            ids(res.uncontrollable_set),
            nothing;
            node = nd,
        )
        plot!(
            p,
            problem;
            plot_region = false,
            observation_region_alpha = 0.0,
            observation_colors = OBSERVATION_COLORS,
            observation_linewidth = 2.0,
        )
        return p
    end
    lims = align_panels!(panels)
    return plot(
        panels...;
        layout = (1, length(panels)),
        size = panel_row_size(lims, length(panels)),
    )
end

save(node_panels(r2.quotient, syn2, "Route 2 ∃"), "spec_route2_synthesis")
save(node_panels(r2.quotient, ver2, "Route 2 ∀"), "spec_route2_verification")

# ── 3-D: the node as a vertical axis ─────────────────────────────────────────────────────────────
# Everything Plots-based must come first — both packages export `plot`.
using CairoMakie

save_mk(mk, name) = begin
    CairoMakie.save(joinpath(@__DIR__, name * ".png"), mk; px_per_unit = 3)
    return println("wrote $name.png")
end

node_z = Dict((1,) => 0.0, (2,) => 1.0)

"Green over red, one layer per node, with the closed loop when there is one."
function augmented(q, res, ttl, zlabel, zt, nz; traj = true)
    mk = CairoMakie.Figure(; size = (900, 700))
    ax = CairoMakie.Axis3(
        mk[1, 1];
        xlabel = "x₁",
        ylabel = "x₂",
        zlabel = zlabel,
        zticks = zt,
        azimuth = 1.2π,
        elevation = 0.16π,
        title = ttl,
    )
    DI.plot_augmented_bisimulation!(
        ax,
        q;
        state_ids = res.uncontrollable_set,
        node_z = nz,
        color_by = LOSING_COLOR,
        alpha = 0.30,
        show_contours = false,
    )
    DI.plot_augmented_bisimulation!(
        ax,
        q;
        state_ids = res.controllable_set,
        node_z = nz,
        color_by = WINNING_COLOR,
        alpha = 0.45,
        show_contours = false,
    )
    traj &&
        res.X !== nothing &&
        DI.plot_augmented_trajectory!(ax, q, res.X, res.M; node_z = nz)
    return mk
end

zt2 = ([0.0, 1.0], ["node (1)", "node (2)"])
save_mk(
    augmented(
        r2.quotient,
        syn2,
        "Route 2 — satisfaction per layer (∃)",
        "memory",
        zt2,
        node_z,
    ),
    "spec_route2_3d_synthesis",
)
save_mk(
    augmented(
        r2.quotient,
        ver2,
        "Route 2 — satisfaction per layer (∀)",
        "memory",
        zt2,
        node_z;
        traj = false,
    ),
    "spec_route2_3d_verification",
)

# Route 1 flat, for the side-by-side: the same certified information with the memory folded away.
zt1 = ([0.0, 1.0], ["", ""])
save_mk(
    augmented(
        r1.quotient,
        syn1,
        "Route 1 — one layer, no memory (∃)",
        "no memory",
        zt1,
        Dict(:clf => 0.0),
    ),
    "spec_route1_3d_synthesis",
)

println(
    """

READ THE TABLE ABOVE FIRST. Route 2 has the larger quotient and the smaller Σ parts; the fixed point
scales with cells, every geometric test scales with parts, and it is the parts that dominate here.
That is the same reason route 2 builds faster, measured on the solve instead of on the build.""",
)
