# ROUTE 2 — build the quotient on the PCLF directly, one partition per node. No determinisation.
#
#     julia --project=test .../A_identical_pieces/demo_route2_direct_on_pclf.jl [dual|primal] [rungs]
#
# The companion of `demo_route1_determinise_first.jl`; run that one first with the SAME arguments.
# Defaults: `dual`, 7 rungs.
#
# Everything below the certificate is identical to route 1's script by construction — same system,
# same two regions, same PCLF, same ladder, same outer level, same tolerance — so the only difference
# between the two is the single line route 1 has and this one does not:
#
#     common = PCLF.build_common_lyapunov(pclf)        # route 1 only
#
# WHAT TO EXPECT, AND WHY THE TWO CONFIGURATIONS DISAGREE
#
#   dual     route 2 builds about 2.4x FEWER cells than route 1 and finishes sooner. Its node commits
#            to the mode it plays next, so it refines under ONE mode instead of two.
#
#   primal   route 2 builds about 1.8x MORE cells and STILL finishes sooner. The graph is complete, so
#            the node constrains nothing about the future, must serve every mode, and does exactly the
#            single node's work — route 2 pays |S| copies of it for nothing. It wins anyway because
#            route 1's cells fragment, and the refinement primitives are superlinear in the number of
#            parts a cell carries.
#
# Same system, same regions, same ladder: only the graph's orientation changes. Route 1 builds nearly
# the same number of cells in both, so the swing is entirely on this side.
#
# WHY GAMMA_X IS HARD-CODED AND TAKEN FROM ROUTE 1. The two arms hold different certificates, so left
# alone they would tile different regions and the comparison would be measuring coverage. The value is
# route 1's own natural level: the baseline keeps the geometry it would have chosen and it is this
# script that is made to conform. Taking it from route 2 instead would flatter route 2.

include(joinpath(dirname(dirname(dirname(@__DIR__))), "common.jl"))

using Printf
import Statistics

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)

const CONFIGS = Dict(
    "dual" =>
        (; dual = true, gamma_x = 1.2747, headline = "expect FEWER cells than route 1"),
    "primal" => (;
        dual = false,
        gamma_x = 1.2456,
        headline = "expect MORE cells than route 1, but simpler",
    ),
)

const WHICH = isempty(ARGS) ? "dual" : ARGS[1]
haskey(CONFIGS, WHICH) ||
    error("first argument must be \"dual\" or \"primal\", got $(repr(WHICH))")
const CFG = CONFIGS[WHICH]
const NB_LEVELS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 7
const GAMMA_X = CFG.gamma_x
const SAMPLES = 3

@printf("configuration: %s De Bruijn k=1, %d rungs\n%s\n\n", WHICH, NB_LEVELS, CFG.headline)

# --- the system: identical to route 1's ----------------------------------------------------------
(; f, X, R1, R2) = two_mode_problem()
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])

# --- the certificate: the SAME one, and the same untimed status ----------------------------------
println("computing the PCLF (the given -- not part of the comparison) ...")
graph = PCLF.generate_DeBruijn_edges(2, 1; dual = CFG.dual)
nodes = sort(collect(graph.verts); by = string)
pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f,
    graph,
    OPTIMIZER;
    Gmats = rotation_templates(nodes; θ = π / 6, mode = :rotation),
    MLF = true,
    verbose = false,
)
@printf("  certified rate %.6f over %d nodes\n", pclf.JSRapprox, length(nodes))
println(
    "  (route 1 certifies the SAME rate — determinisation preserves it, which is what makes a",
)
println("   cost comparison between the two meaningful)\n")

# --- ROUTE 2 -------------------------------------------------------------------------------------
# The PCLF goes straight in. Each node keeps its own piece and refines only under the modes its own
# outgoing edges carry.
route2() = build_quotient(
    problem,
    pclf;
    atol = 1e-3,
    nb_levels = NB_LEVELS,
    ΓX = GAMMA_X,
    max_slices = NB_LEVELS,
    print_level = 0,
).quotient

println(
    "warming up (the first call compiles the whole pipeline; timing it would measure Julia) ...",
)
route2()

println("timing route 2 over $SAMPLES runs ...")
times = Float64[]
for _ in 1:SAMPLES
    GC.gc()
    t0 = time_ns()
    route2()
    push!(times, (time_ns() - t0) / 1e9)
end

q = route2()
parts, faces = PCQ.cell_complexities(q)
st = PCQ.bisimulation_stats(q)

println("\n", "="^76)
@printf("ROUTE 2 — straight onto the PCLF   [%s, %d rungs]\n", WHICH, NB_LEVELS)
println("="^76)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)\n\n",
    "BUILD TIME",
    Statistics.median(times),
    minimum(times),
    maximum(times)
)
@printf("%-30s %d\n", "cells", st[:num_states])
@printf("%-30s %d\n", "transitions", st[:num_transitions])
@printf("%-30s %d\n", "slices", st[:num_slices])
@printf("%-30s %d\n", "nodes carrying a partition", st[:num_nodes])
@printf("%-30s %d\n", "MAX facets in one cell", maximum(faces))
@printf("%-30s %d\n", "MAX parts in one cell", maximum(parts))

# --- what this route covers, node by node and in total -------------------------------------------
# Each node tiles its OWN piece's sublevel set, and those overlap in the plane, so the per-node areas
# do not add up to anything meaningful -- the covered region is their UNION, reported last.
println()
node_sets = LazySets.HPolytope[]
tiled_areas = Float64[]          # a vector, not an accumulator: a top-level `for` opens a new scope
for nd in quotient_nodes(q)
    ids = [k for (k, v) in q.states if v.node == nd]
    a = PCQ.get_volume(q, ids; backend = CDDLib.Library())
    push!(tiled_areas, a)
    @printf("%-30s %.4f\n", "  area tiled by node $nd", a)
    S = PCLF.get_sublevel_set(pclf.pieces[nd], GAMMA_X)
    append!(node_sets, S isa LazySets.UnionSetArray ? S.array : [UT._as_hpolytope(S)])
end
# COMPARE LIKE WITH LIKE. The first line is the region the certificate defines -- the union of the
# node pieces' sublevel sets -- and it is what route 1's "region the certificate defines" must be set
# against. The second is what was actually tiled, summed over nodes; it is larger than the union
# because the nodes OVERLAP in the plane, and smaller per node than the ideal because `atol` erodes
# every cut. Neither of those two effects is a difference of domain, and mixing the two lines across
# routes invents one.
@printf(
    "%-30s %.4f  <- compare with route 1's SAME line\n",
    "region the certificate defines",
    UT.get_volume(UT.semilinear_set(node_sets); backend = CDDLib.Library())
)
@printf(
    "%-30s %.4f  (nodes overlap, so this is not an area)\n",
    "sum of tiled areas",
    sum(tiled_areas)
)

# --- the picture, drawn strictly AFTER every measurement -----------------------------------------
# The nodes are drawn SEPARATELY: merged they would superimpose |S| partitions over one plane and show
# neither, and separated they ARE the argument — route 1's single plane against these.
gr()
figdir = @__DIR__
mkpath(figdir)
panels = map(quotient_nodes(q)) do nd
    n_here = count(s -> s.node == nd, values(q.states))
    p = plot(;
        aspect_ratio = :equal,
        legend = false,
        title = "route 2 — node $nd\n$n_here cells",
        titlefontsize = 9,
    )
    plot!(
        p,
        q;
        what = :states,
        by = :state,
        node = nd,
        show_contours = true,
        linewidth = 0.3,
        fillalpha = 0.9,
        merge_series = false,
    )
    return p
end
lims = align_panels!(panels)
fig = plot(
    panels...;
    layout = (1, length(panels)),
    size = panel_row_size(lims, length(panels)),
    plot_title = "route 2 — $(st[:num_states]) cells over $(length(panels)) nodes, worst cell $(maximum(parts)) parts",
    plot_titlefontsize = 10,
)
out = joinpath(figdir, "demo_route2_$(WHICH).png")
savefig(fig, out)
println("\nwrote $out")

println(
    """
COMPARE WITH ROUTE 1, on both counts.

  cells          on `dual` this arm has fewer, on `primal` more. The count alone does not decide.
  parts per cell this arm's worst cell is a handful of convex parts; route 1's runs into the hundreds
                 on the primal graph, because its induced common there is a non-convex UNION and every
                 cell inherits it. Every set difference, pre-image and emptiness test pays for those.
  covered area   on `primal` the two agree (route 1's common IS the union). On `dual` route 1's common
                 is the INTERSECTION and it certifies strictly LESS than the union reported above — a
                 point lying in one node's piece but not the other is certified here and not there.

One honest caveat: a ratio from one run of each on a busy machine is indicative, not a measurement.
The campaign's figures come from interleaved rounds with the ratio taken pairwise per round — see
plan.md §5, and `1_comparison.jl` for the full matrix.""",
)
