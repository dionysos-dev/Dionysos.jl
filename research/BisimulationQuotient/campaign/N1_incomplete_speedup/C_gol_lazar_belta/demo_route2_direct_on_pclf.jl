# GOL–LAZAR–BELTA, ROUTE 2 — build the quotient straight on the PCLF, then synthesise AND verify.
#
#     julia --project=test .../C_gol_lazar_belta/demo_route2_direct_on_pclf.jl [rungs]
#
# The companion of `demo_route1_determinise_first.jl`; run that one first with the SAME number of rungs.
#
# Everything below the certificate is identical to route 1's script by construction — their dynamics,
# their working set, their three regions, their formula, their initial point, the same PCLF, the same
# ladder, the same outer level, the same tolerance. The only difference is the single line route 1 has
# and this one does not:
#
#     common = PCLF.build_common_lyapunov(pclf)        # route 1 only
#
# WHAT TO EXPECT
#
# About 1.6x MORE cells than route 1, and several times less wall time. The count and the clock point
# in opposite directions here, which is the clearest statement the campaign can make: **cell count is
# the wrong cost proxy**. Route 1's induced common on this complete graph is a non-convex UNION of 17
# parts; every one of its cells inherits that fragmentation, and every set difference, pre-image and
# emptiness test pays for the parts. This arm keeps two simple pieces instead.
#
# BOTH QUANTIFIERS ARE SOLVED, because the paper's example poses both and because a cheaper quotient
# is only worth having if it answers both questions the same way:
#
#   synthesis (∃)     the modes are the CONTROLLER's — does some switching signal satisfy φ?
#   verification (∀)  the modes are the ENVIRONMENT's — does every signal satisfy it?
#
# Compare the two green sets against route 1's. If they disagree, every cost number here is void, so
# those figures are the check rather than the decoration.

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
const NB_LEVELS = isempty(ARGS) ? 5 : parse(Int, ARGS[1])
const GAMMA_X = 5.5202      # route 1's outer level, so both arms tile the same set
const SAMPLES = 3

@printf(
    "Gol–Lazar–Belta Example 3.1 — route 2, primal De Bruijn k=1, %d rungs\n\n",
    NB_LEVELS
)

# --- their problem: identical to route 1's --------------------------------------------------------
(; f, X, R1, R2, R3) = gol_lazar_belta_problem()
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])
φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)     # point `a` of the paper

# --- the certificate: the SAME one, and the same untimed status -----------------------------------
println("computing the PCLF (the given — not part of the comparison) ...")
graph = PCLF.generate_DeBruijn_edges(2, 1; dual = false)
nodes = sort(collect(graph.verts); by = string)
pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    graph,
    OPTIMIZER,
    PCLF.conic_partitions_dict_2d(2, nodes);
    MLF = true,
)
@printf("  certified rate %.6f over %d nodes\n\n", pclf.JSRapprox, length(nodes))

# --- ROUTE 2: the PCLF goes straight in -----------------------------------------------------------
build() = build_quotient(
    problem,
    pclf;
    atol = 1e-3,
    nb_levels = NB_LEVELS,
    ΓX = GAMMA_X,
    max_slices = NB_LEVELS,
    print_level = 0,
)

println(
    "warming up (the first call compiles the pipeline; timing it would measure Julia) ...",
)
build()

println("timing the build over $SAMPLES runs ...")
t_build = Float64[]
for _ in 1:SAMPLES
    GC.gc()
    t0 = time_ns()
    build()
    push!(t_build, (time_ns() - t0) / 1e9)
end

built = build()
q, D = built.quotient, built.D
parts, faces = PCQ.cell_complexities(q)
st = PCQ.bisimulation_stats(q)

# --- the specification, both quantifiers, timed separately ----------------------------------------
regions = Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3)
ap_to_obs = Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3)

"Solve φ on the quotient. `system` decides who owns the modes, and therefore the quantifier."
solver(system) =
    () -> synthesize_cosafe_ltl(
        system,
        q,
        Dionysos.spot_stepper(φ),
        regions,
        ap_to_obs,
        x0;
        print_level = 0,
    )

solve = solver(f)                                                        # ∃: the controller switches
verify = solver(ST.with_switching(f, HybridSystems.AutonomousSwitching()))  # ∀: the environment does

function timed(run, label)
    println("warming up the $label solve ...")
    run()
    println("timing the $label solve over $SAMPLES runs ...")
    ts = Float64[]
    for _ in 1:SAMPLES
        GC.gc()
        t0 = time_ns()
        run()
        push!(ts, (time_ns() - t0) / 1e9)
    end
    return run(), ts
end

result, t_solve = timed(solve, "synthesis")
verification, t_verify = timed(verify, "verification")

med(v) = Statistics.median(v)
println("\n", "="^76)
@printf("ROUTE 2 — straight onto the PCLF   [Gol–Lazar–Belta, %d rungs]\n", NB_LEVELS)
println("="^76)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)\n",
    "BUILD TIME",
    med(t_build),
    minimum(t_build),
    maximum(t_build)
)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)\n",
    "co-safe synthesis (∃)",
    med(t_solve),
    minimum(t_solve),
    maximum(t_solve)
)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)\n\n",
    "verification (∀)",
    med(t_verify),
    minimum(t_verify),
    maximum(t_verify)
)
@printf("%-30s %d\n", "cells", st[:num_states])
@printf("%-30s %d\n", "transitions", st[:num_transitions])
@printf("%-30s %d\n", "slices", st[:num_slices])
@printf("%-30s %d\n", "nodes carrying a partition", st[:num_nodes])
@printf("%-30s %d\n", "MAX facets in one cell", maximum(faces))
@printf("%-30s %d\n", "MAX parts in one cell", maximum(parts))
@printf("%-30s %d\n", "cells won (synthesis, ∃)", length(result.controllable_set))
@printf("%-30s %d\n", "cells verified (∀)", length(verification.controllable_set))

# --- the pictures, drawn strictly AFTER every measurement ----------------------------------------
figdir = @__DIR__
mkpath(figdir)

# The nodes are drawn separately: merged they would superimpose |S| partitions over one plane and
# show neither, and separated they are the argument.
panels = map(quotient_nodes(q)) do nd
    n_here = count(s -> s.node == nd, values(q.states))
    p = plot(;
        aspect_ratio = :equal,
        legend = false,
        title = "node $nd — $n_here cells",
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
savefig(
    plot(
        panels...;
        layout = (1, length(panels)),
        size = panel_row_size(lims, length(panels)),
        plot_title = "GLB route 2 — $(st[:num_states]) cells over $(length(panels)) nodes, worst cell $(maximum(parts)) parts",
        plot_titlefontsize = 10,
    ),
    joinpath(figdir, "demo_route2_quotient.png"),
)
savefig(
    plot_synthesis_result(q, result, problem, "GLB route 2 — SYNTHESIS (∃) of φ"),
    joinpath(figdir, "demo_route2_synthesis.png"),
)
savefig(
    plot_synthesis_result(q, verification, problem, "GLB route 2 — VERIFICATION (∀) of φ"),
    joinpath(figdir, "demo_route2_verification.png"),
)
println("\nwrote demo_route2_{quotient,synthesis,verification}.png in $figdir")

println(
    """
COMPARE WITH ROUTE 1 on four things.

  time           this arm should be several times faster to build.
  cells          it has MORE of them. The count alone decides nothing.
  parts per cell its worst cell holds a handful of convex parts; route 1's holds many more, because
                 the induced common it works on is a non-convex union. That is where the time goes.
  green sets     BOTH quantifiers, against route 1's. Same formula, same answer, two very different
                 partitions — if they disagree, the cost numbers above are void.

One caveat: a ratio from one run of each on a busy machine is indicative, not a measurement. The
campaign's figures come from interleaved rounds with the ratio taken pairwise per round (plan.md §5).""",
)
