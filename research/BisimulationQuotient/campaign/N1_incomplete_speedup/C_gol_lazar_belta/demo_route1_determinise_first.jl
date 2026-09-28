# GOL–LAZAR–BELTA, ROUTE 1 — determinise the PCLF, then build the quotient and synthesise on it.
#
#     julia --project=test .../C_gol_lazar_belta/demo_route1_determinise_first.jl [rungs]
#     julia --project=test .../C_gol_lazar_belta/demo_route2_direct_on_pclf.jl [rungs]
#
# Default 5 rungs. Use the SAME value for both, or the comparison is meaningless.
#
# THE PROBLEM IS THEIRS, THE CERTIFICATE IS OURS
#
# Dynamics, working set, the three observation regions, the co-safe formula and the initial point are
# Example 3.1 of Gol, Ding, Lazar & Belta (arXiv:1208.5471) — a benchmark a reader already knows.
#
#     φ = (!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D)
#
# What is NOT compared against is their published certificate, a common Lyapunov function
# (‖Lx‖_∞, ρ = 0.94). A certificate determines its own slice family and terminal set, so two
# different certificates tile two different regions and their cell counts answer different
# questions. The comparison here is the campaign's — ONE path-complete certificate from an oracle,
# then route 1 against route 2 on it — posed on their problem.
#
# WHAT TO EXPECT (primal De Bruijn k=1, certified rate 0.867388)
#
# Route 2 builds about 1.6x MORE cells than this script and finishes 7x to 9x sooner. That is the
# whole point: the induced common on a complete graph is a non-convex UNION — 17 parts here — every
# cell of the quotient inherits that fragmentation, and the refinement primitives are superlinear in
# the number of parts a cell carries. Cell count is the wrong cost proxy.
#
# WHY GAMMA_X IS HARD-CODED. Each arm would otherwise derive its outer level from its own certificate,
# start from a different level, tile a different region, and their counts would measure coverage
# instead of efficiency. The value is route 1's own natural level, so the baseline keeps the geometry
# it would have chosen and route 2 is made to conform; taking it from route 2 would flatter route 2.

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
const GAMMA_X = 5.5202      # route 1's own outer level on this problem, at any depth
const SAMPLES = 3

@printf(
    "Gol–Lazar–Belta Example 3.1 — route 1, primal De Bruijn k=1, %d rungs\n\n",
    NB_LEVELS
)

# --- their problem -------------------------------------------------------------------------------
(; f, X, R1, R2, R3) = gol_lazar_belta_problem()
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])
φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)     # point `a` of the paper

# --- our certificate: GIVEN, and deliberately not timed -------------------------------------------
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

quotient_of(cert) = build_quotient(
    problem,
    cert;
    atol = 1e-3,
    nb_levels = NB_LEVELS,
    ΓX = GAMMA_X,
    max_slices = NB_LEVELS,
    print_level = 0,
)

println(
    "warming up (the first call compiles the pipeline; timing it would measure Julia) ...",
)
quotient_of(PCLF.build_common_lyapunov(pclf))

# The two phases are timed separately, because "determinising is slow" and "what determinising
# returns is slow to build on" are different claims and only the second is true.
println("timing the build over $SAMPLES runs, phase by phase ...")
t_det, t_build = Float64[], Float64[]
for _ in 1:SAMPLES
    GC.gc()
    t0 = time_ns()
    common = PCLF.build_common_lyapunov(pclf)
    t1 = time_ns()
    quotient_of(common)
    t2 = time_ns()
    push!(t_det, (t1 - t0) / 1e9)
    push!(t_build, (t2 - t1) / 1e9)
end

common_cert = PCLF.build_common_lyapunov(pclf)
common_parts = let
    S = PCLF.get_sublevel_set(common_cert.pieces[:clf], 1.0; atol = 1e-6)
    S isa LazySets.UnionSetArray ? length(S.array) : 1
end
built = quotient_of(common_cert)
q, D = built.quotient, built.D
parts, faces = PCQ.cell_complexities(q)
st = PCQ.bisimulation_stats(q)

# --- the specification, timed separately ---------------------------------------------------------
regions = Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3)
ap_to_obs = Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3)
solve() = synthesize_cosafe_ltl(
    f,
    q,
    Dionysos.spot_stepper(φ),
    regions,
    ap_to_obs,
    x0;
    print_level = 0,
)
println("warming up the co-safe solve ...")
solve()
println("timing the co-safe solve over $SAMPLES runs ...")
t_solve = Float64[]
for _ in 1:SAMPLES
    GC.gc()
    t0 = time_ns()
    solve()
    push!(t_solve, (time_ns() - t0) / 1e9)
end
result = solve()

# --- the SAME formula under VERIFICATION (∀) ------------------------------------------------------
# Synthesis asks whether SOME switching signal satisfies φ — the modes are the controller's, and it
# picks. Verification asks whether EVERY signal does: the modes become the environment's, and the
# adversary picks. Same quotient, same formula, opposite quantifier, and a much smaller green set.
# Both are reported because the paper's example poses both, and because a quotient that is cheaper to
# build is only worth having if it answers both questions identically.
f_forall = ST.with_switching(f, HybridSystems.AutonomousSwitching())
verify() = synthesize_cosafe_ltl(
    f_forall,
    q,
    Dionysos.spot_stepper(φ),
    regions,
    ap_to_obs,
    x0;
    print_level = 0,
)
println("warming up the verification solve ...")
verify()
println("timing the verification solve over $SAMPLES runs ...")
t_verify = Float64[]
for _ in 1:SAMPLES
    GC.gc()
    t0 = time_ns()
    verify()
    push!(t_verify, (time_ns() - t0) / 1e9)
end
verification = verify()

med(v) = Statistics.median(v)
tot = t_det .+ t_build
println("\n", "="^76)
@printf("ROUTE 1 — determinise first   [Gol–Lazar–Belta, %d rungs]\n", NB_LEVELS)
println("="^76)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)  %5.1f %%\n",
    "1. determinise the PCLF",
    med(t_det),
    minimum(t_det),
    maximum(t_det),
    100 * med(t_det) / med(tot)
)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)  %5.1f %%\n",
    "2. build the quotient on it",
    med(t_build),
    minimum(t_build),
    maximum(t_build),
    100 * med(t_build) / med(tot)
)
@printf("%-30s %8.3f s\n", "BUILD TOTAL", med(tot))
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)\n\n",
    "3. co-safe synthesis (∃)",
    med(t_solve),
    minimum(t_solve),
    maximum(t_solve)
)
@printf(
    "%-30s %8.3f s  (min %.3f, max %.3f)\n\n",
    "4. verification (∀)",
    med(t_verify),
    minimum(t_verify),
    maximum(t_verify)
)
@printf("%-30s %d\n", "parts in the induced common", common_parts)
@printf("%-30s %d\n", "cells", st[:num_states])
@printf("%-30s %d\n", "transitions", st[:num_transitions])
@printf("%-30s %d\n", "slices", st[:num_slices])
@printf("%-30s %d\n", "MAX facets in one cell", maximum(faces))
@printf("%-30s %d\n", "MAX parts in one cell", maximum(parts))
@printf("%-30s %d\n", "cells won (synthesis, ∃)", length(result.controllable_set))
@printf("%-30s %d\n", "cells verified (∀)", length(verification.controllable_set))

# --- the pictures, drawn strictly AFTER every measurement ----------------------------------------
figdir = @__DIR__
mkpath(figdir)

fig = plot(;
    aspect_ratio = :equal,
    legend = false,
    title = "GLB route 1 — the induced common\n$(st[:num_states]) cells on ONE plane, worst cell $(maximum(parts)) parts",
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
savefig(fig, joinpath(figdir, "demo_route1_quotient.png"))

savefig(
    plot_synthesis_result(q, result, problem, "GLB route 1 — SYNTHESIS (∃) of φ"),
    joinpath(figdir, "demo_route1_synthesis.png"),
)
savefig(
    plot_synthesis_result(q, verification, problem, "GLB route 1 — VERIFICATION (∀) of φ"),
    joinpath(figdir, "demo_route1_verification.png"),
)
println("\nwrote demo_route1_{quotient,synthesis,verification}.png in $figdir")

println(
    """
READ THE SPLIT FIRST: determinising costs essentially nothing. What costs is building on what it
returned — an induced common that is a $(common_parts)-part non-convex union on this complete graph.

Now run `demo_route2_direct_on_pclf.jl $NB_LEVELS`. It builds MORE cells than this and finishes several times
sooner, and it certifies the same formula.""",
)
