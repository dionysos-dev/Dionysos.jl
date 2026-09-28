# ROUTE 1 — determinise the PCLF into a common Lyapunov function, then build the quotient on it.
#
#     julia --project=test .../A_identical_pieces/demo_route1_determinise_first.jl [dual|primal] [rungs]
#     julia --project=test .../A_identical_pieces/demo_route2_direct_on_pclf.jl    [dual|primal] [rungs]
#
# Defaults: `dual`, 7 rungs. Use the SAME arguments for both scripts or the comparison is meaningless.
#
# TWO CONFIGURATIONS, TWO DIFFERENT MECHANISMS, ONE SYSTEM
#
# Everything is held fixed except the orientation of the De Bruijn graph -- same modes, same working
# set, same two observation regions, same ladder, same template. Route 1 builds nearly the SAME number
# of cells either way (about 3 500 and 3 700); it is route 2 that swings, and in opposite directions:
#
#   dual     the node commits to the mode it will play NEXT, so it refines under ONE mode.
#            Route 2 builds ~2.4x FEWER cells and is faster.
#
#   primal   the node records the mode just PLAYED, which constrains nothing about the future, so it
#            must still serve every mode and does exactly the single node's work. Route 2 pays |S|
#            copies of that and builds ~1.8x MORE cells -- and is still faster, because route 1's
#            cells fragment: the induced common on a complete graph is a non-convex UNION of the node
#            pieces and every cell inherits it.
#
# So "route 2 builds a smaller quotient" is true on one graph and false on the other, while "route 2
# finishes sooner" is true on both. That is the point worth taking away, and it is why cell count is
# the wrong cost proxy.
#
# WHAT THE TWO SCRIPTS SHARE, SO THAT THE COMPARISON IS FAIR
#
#   same system        `two_mode_problem`, two observation regions
#   same certificate   ONE path-complete Lyapunov function, computed identically in both and NOT
#                      timed: it is the given, an oracle's output, not part of what is compared
#   same ladder        the same number of rungs and the same outer level GAMMA_X
#   same tolerance     atol = 1e-3
#
# WHY GAMMA_X IS HARD-CODED. The construction normally derives its outer level from the certificate it
# is handed, and the two arms hold different certificates -- route 1 the induced common, route 2 the
# PCLF -- so left alone they would start from different levels, tile different regions, and their cell
# counts would measure coverage rather than efficiency. The values below are route 1's own natural
# levels, so the baseline keeps the geometry it would have chosen and route 2 is made to conform.
# Choosing them the other way round would flatter route 2.

include(joinpath(dirname(dirname(dirname(@__DIR__))), "common.jl"))

using Printf
import Statistics

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)

const CONFIGS = Dict(
    "dual" => (;
        dual = true,
        gamma_x = 1.2747,
        headline = "route 2 builds FEWER cells (about 2.4x fewer)",
    ),
    "primal" => (;
        dual = false,
        gamma_x = 1.2456,
        headline = "route 2 builds MORE cells (about 1.8x more) but far SIMPLER ones",
    ),
)

const WHICH = isempty(ARGS) ? "dual" : ARGS[1]
haskey(CONFIGS, WHICH) ||
    error("first argument must be \"dual\" or \"primal\", got $(repr(WHICH))")
const CFG = CONFIGS[WHICH]

# Each rung multiplies the partition, so this is the "make the quotient bigger" knob. At 7 rungs the
# terminal set is about a tenth of the outer level, which leaves the ladder room to be visible; at 5
# it is a quarter and the innermost cell dominates the picture. Above 8, route 1 runs into minutes.
const NB_LEVELS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 7
const GAMMA_X = CFG.gamma_x
const SAMPLES = 3

@printf("configuration: %s De Bruijn k=1, %d rungs\n%s\n\n", WHICH, NB_LEVELS, CFG.headline)

# --- the system ----------------------------------------------------------------------------------
(; f, X, R1, R2) = two_mode_problem()
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])

# --- the certificate: GIVEN, and deliberately not timed ------------------------------------------
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
@printf("  certified rate %.6f over %d nodes\n\n", pclf.JSRapprox, length(nodes))

quotient_from(cert) = build_quotient(
    problem,
    cert;
    atol = 1e-3,
    nb_levels = NB_LEVELS,
    ΓX = GAMMA_X,
    max_slices = NB_LEVELS,
    print_level = 0,
).quotient

println(
    "warming up (the first call compiles the whole pipeline; timing it would measure Julia) ...",
)
quotient_from(PCLF.build_common_lyapunov(pclf))

# The two phases are timed SEPARATELY, because "determinising is expensive" and "what determinising
# produces is expensive to build on" are different claims and only the second is true. The subset
# construction is combinatorics on a two-node graph and costs essentially nothing.
println("timing route 1 over $SAMPLES runs, phase by phase ...")
t_determinise, t_build = Float64[], Float64[]
for _ in 1:SAMPLES
    GC.gc()
    t0 = time_ns()
    common = PCLF.build_common_lyapunov(pclf)
    t1 = time_ns()
    quotient_from(common)
    t2 = time_ns()
    push!(t_determinise, (t1 - t0) / 1e9)
    push!(t_build, (t2 - t1) / 1e9)
end
times = t_determinise .+ t_build

common_cert = PCLF.build_common_lyapunov(pclf)
common_parts = let
    S = PCLF.get_sublevel_set(common_cert.pieces[:clf], 1.0; atol = 1e-6)
    S isa LazySets.UnionSetArray ? length(S.array) : 1
end
q = quotient_from(common_cert)
parts, faces = PCQ.cell_complexities(q)
st = PCQ.bisimulation_stats(q)
med(v) = Statistics.median(v)

println("\n", "="^76)
@printf("ROUTE 1 — determinise first, then build   [%s, %d rungs]\n", WHICH, NB_LEVELS)
println("="^76)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)  %5.1f %%\n",
    "1. determinise the PCLF",
    med(t_determinise),
    minimum(t_determinise),
    maximum(t_determinise),
    100 * med(t_determinise) / med(times)
)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)  %5.1f %%\n",
    "2. build the quotient on it",
    med(t_build),
    minimum(t_build),
    maximum(t_build),
    100 * med(t_build) / med(times)
)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)\n\n",
    "TOTAL",
    med(times),
    minimum(times),
    maximum(times)
)
@printf("%-30s %d\n", "parts in the induced common", common_parts)
@printf("%-30s %d\n", "cells", st[:num_states])
@printf("%-30s %d\n", "transitions", st[:num_transitions])
@printf("%-30s %d\n", "slices", st[:num_slices])
@printf("%-30s %d\n", "MAX facets in one cell", maximum(faces))
@printf("%-30s %d\n", "MAX parts in one cell", maximum(parts))
# TWO areas, and conflating them is a mistake this script made once. The region the CERTIFICATE
# defines is what the two routes should be compared on -- it is the theory's claim. The area actually
# TILED is smaller, because `atol` insets every cut and the arm making more cuts loses more; that gap
# is a cost of fine cutting, not a difference of domain. Comparing one route's tiled area against the
# other's certificate region mixes the two and manufactures a difference out of erosion.
region = PCLF.get_sublevel_set(common_cert.pieces[:clf], GAMMA_X; atol = 1e-6)
@printf(
    "%-30s %.4f\n",
    "region the certificate defines",
    UT.get_volume(region; backend = CDDLib.Library())
)
@printf(
    "%-30s %.4f\n",
    "area actually tiled",
    PCQ.get_volume(q, keys(q.states); backend = CDDLib.Library())
)

# --- the picture, drawn strictly AFTER every measurement -----------------------------------------
gr()
figdir = @__DIR__
mkpath(figdir)
fig = plot(;
    aspect_ratio = :equal,
    legend = false,
    title = "route 1 — the induced common\n$(st[:num_states]) cells on ONE plane, worst cell $(maximum(parts)) parts",
    titlefontsize = 9,
)
plot!(
    fig,
    q;
    what = :states,
    by = :state,
    show_contours = true,
    linewidth = 0.3,
    fillalpha = 0.9,
    merge_series = false,
)
lims = align_panels!([fig])
plot!(fig; size = panel_row_size(lims, 1))
out = joinpath(figdir, "demo_route1_$(WHICH).png")
savefig(fig, out)
println("\nwrote $out")

println(
    """
READ THE SPLIT FIRST. Determinising is essentially free — it is a subset construction on a two-node
graph with no geometry in it. Practically the whole bill is phase 2, building the quotient on what
that step produced. The claim is therefore NOT "determinising is slow"; it is that what determinisation
RETURNS is expensive to build on.

ON COVERAGE, and it depends on the graph:
  primal (complete)     the induced common is min_i V_i, whose sublevel set is the UNION of the node
                        pieces — exactly what route 2's nodes tile between them. SAME region.
  dual (co-complete)    the induced common is max_i V_i, the INTERSECTION — strictly smaller. Route 1
                        then certifies LESS than route 2, not merely more slowly.
Compare the two "region the certificate defines" lines -- never a tiled area against a certificate
region, which is how this script once manufactured a 14 % difference out of `atol` erosion.

Now run `demo_route2_direct_on_pclf.jl $WHICH $NB_LEVELS`.""",
)
