# GLB, the POSTER's configuration exactly — route 1 against route 2, fully characterised.
#
#     julia --project=test .../C_gol_lazar_belta/poster_configuration.jl        SAMPLES=n
#
# `examples/gol_lazar_belta_pclf.jl` produced the poster figures. This reproduces its settings
# verbatim and reports, for BOTH routes, everything needed to decide whether the example can carry
# the paper: the certificate, the graph, build time, and the full per-cell complexity distribution.
#
# THE POSTER'S SETTINGS, AND WHERE THE CAMPAIGN RUNS DIFFERED
#
#                        poster (here)                  campaign runs
#   working set X        the ±5.9 box                   their {‖Lx‖∞ ≤ 10}
#   ladder               natural stop, max_slices = 50  nb_levels forced, ΓX pinned
#   atol                 1e-4                           1e-3
#   graph, template      identical — primal De Bruijn order 1, conic order 2, MLF
#
# `conic_partitions_dict_2d(order, nodes)` is literally `Dict(id => conic_partitions_2d(order))`, so
# the certificate is the same object in both; only the domain, the ladder and the tolerance moved.
#
# ONE THING TO READ BEFORE TRUSTING THE COMPARISON. The poster lets **each arm derive its own** outer
# level from its own certificate, and the induced common's gauge is not the PCLF pieces' gauge — so
# the two arms may tile different regions, and a cost ratio between them would then be measuring
# coverage. The script prints each arm's τX and slice count. If they agree the comparison is
# like-for-like as it stands; if they do not, pin ΓX before quoting a ratio.

include(joinpath(dirname(dirname(dirname(@__DIR__))), "common.jl"))

using Printf
import Statistics

const SAMPLES = parse(Int, get(ENV, "SAMPLES", "1"))
const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)

# ── the problem, exactly as the poster takes it ──────────────────────────────────────────────────
(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
regions = [UT._as_hpolytope(R) for R in (R1, R2, R3)]

# ── the certificate, exactly as the poster builds it ─────────────────────────────────────────────
graph = PCLF.generate_DeBruijn_edges(2, 1)          # primal De Bruijn, order 1
nodes = sort(collect(graph.verts); by = string)
partition = PCLF.conic_partitions_2d(2)             # conic partition of the plane, order 2
pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    graph,
    OPTIMIZER,
    Dict(nd => partition for nd in nodes);
    MLF = true,
)
common = PCLF.build_common_lyapunov(pclf)

println("="^84)
println("THE CERTIFICATE")
println("="^84)
@printf(
    "graph            primal De Bruijn, order 1 — %d nodes, %d edges\n",
    length(graph.verts),
    length(graph.edges)
)
for (u, v, m) in sort(graph.edges; by = e -> (string(e[1]), e[3]))
    @printf("                 %s --mode %d--> %s\n", string(u), Int(m), string(v))
end
@printf(
    "                 complete: %s   co-complete: %s\n",
    PCLF.is_complete(graph, 1:2),
    PCLF.is_co_complete(graph, 1:2)
)
@printf(
    "template         conic partition of the plane, order 2 — %d cones, shared by both nodes\n",
    length(partition)
)
@printf("certified rate   %.6f\n", pclf.JSRapprox)
for nd in nodes
    p = pclf.pieces[nd]
    @printf("  piece %-6s  %d template rows\n", string(nd), size(p.G, 1))
end
let S = PCLF.get_sublevel_set(common.pieces[:clf], 1.0; atol = 1e-6)
    n = S isa LazySets.UnionSetArray ? length(S.array) : 1
    @printf(
        "induced common   %d convex part(s) at γ = 1 — this is what route 1 builds on\n",
        n
    )
end

# ── the two arms, each with the poster's own natural stop ────────────────────────────────────────
build(cert; ΓX = nothing, nb_levels = nothing) = build_quotient(
    problem,
    cert;
    atol = 1e-4,
    level_tol = 1e-2,
    max_slices = nb_levels === nothing ? 50 : nb_levels,
    ΓX = ΓX,
    nb_levels = nb_levels,
    print_level = 0,
)

"Every per-cell figure the paper might quote, not just the maxima."
function stats(q)
    parts, faces = PCQ.cell_complexities(q)
    st = PCQ.bisimulation_stats(q)
    return (;
        cells = st[:num_states],
        transitions = st[:num_transitions],
        slices = st[:num_slices],
        nodes = st[:num_nodes],
        sum_p = sum(parts),
        mean_p = Statistics.mean(parts),
        med_p = Statistics.median(parts),
        max_p = maximum(parts),
        sum_f = sum(faces),
        mean_f = Statistics.mean(faces),
        med_f = Statistics.median(faces),
        max_f = maximum(faces),
    )
end

function arm(name, cert; ΓX = nothing, nb_levels = nothing)
    println("\n", "-"^84)
    println(name)
    println("-"^84)
    print("  warming up ... ")
    b = build(cert; ΓX = ΓX, nb_levels = nb_levels)
    println("done")
    ts = Float64[]
    for k in 1:SAMPLES
        GC.gc()
        t0 = time_ns()
        b = build(cert; ΓX = ΓX, nb_levels = nb_levels)
        push!(ts, (time_ns() - t0) / 1e9)
        @printf("  round %d: %.3f s\n", k, ts[k])
    end
    s = stats(b.quotient)
    τ = MOI.get(b.optimizer, MOI.RawOptimizerAttribute("Γ"))
    @printf("  τX = %.4f   τD = %.4f   %d slices\n", maximum(τ), minimum(τ), s.slices)

    # The natural stop is *defined* by this, but the set it returned is what the labelling uses.
    Dparts = b.D isa LazySets.UnionSetArray ? b.D.array : [UT._as_hpolytope(b.D)]
    ok = all(LazySets.isdisjoint(P, R) for P in Dparts, R in regions)
    @printf("  D: %d part(s), clears R1/R2/R3: %s\n", length(Dparts), ok)
    return (; name, t = ts, s, τ)
end

println("\n\n", "#"^84)
println("PASS 1 — the poster's own settings: each arm derives its OWN outer level")
println("#"^84)
r1 = arm("ROUTE 1 — determinise first, then build on the induced common", common)
r2 = arm("ROUTE 2 — build directly on the PCLF, one partition per node", pclf)

med(v) = Statistics.median(v)

# ── PASS 2, only if pass 1 showed the arms starting from different levels ────────────────────────
# Route 1's own natural level and slice count are imposed on route 2: the baseline keeps the geometry
# it would have chosen, and the challenger is made to conform. Taking them from route 2 would flatter
# route 2. Without this the cell counts measure how much of the plane each arm covered.
const PIN_Γ = maximum(r1.τ)
const PIN_N = r1.s.slices
p1, p2 = nothing, nothing
if !isapprox(maximum(r1.τ), maximum(r2.τ); rtol = 0.02)
    println("\n\n", "#"^84)
    @printf(
        "PASS 2 — same ladder for both: ΓX = %.4f, %d slices (route 1's own)\n",
        PIN_Γ,
        PIN_N
    )
    println("#"^84)
    p1 = arm("ROUTE 1 — pinned ladder", common; ΓX = PIN_Γ, nb_levels = PIN_N)
    p2 = arm("ROUTE 2 — pinned ladder", pclf; ΓX = PIN_Γ, nb_levels = PIN_N)
end

# ── PASS 3 (opt-in: PASS3=1) — the level at which EVERY node covers the working set ──────────────
# Pass 2 pins at the level where the UNION of the pieces covers X, which is what the induced common
# needs. On a complete graph that union IS route 1's sublevel set, so both arms tile the same region
# at any common Γ — but an individual node's piece need not cover X there, and the universal
# semantics quantifies over every node. `compute_tau_X(pclf, ·)` takes the max over pieces, which is
# exactly the level that fixes it.
q1, q2 = nothing, nothing
if get(ENV, "PASS3", "0") == "1"
    Xh = UT._as_hpolytope(X)
    Γ3 = PCQ.compute_tau_X(pclf, Xh)
    γ = pclf.JSRapprox
    N3 = something(
        findfirst(
            k ->
                PCQ.all_nodes_clear_regions(pclf, Γ3 * γ^(k - 1), regions; tol = 1e-2) && PCQ.all_nodes_clear_regions(
                    common,
                    Γ3 * γ^(k - 1),
                    regions;
                    tol = 1e-2,
                ),
            1:60,
        ),
        -1,
    )
    N3 > 0 || error("no region-free terminal level within 60 rungs at ΓX = $Γ3")
    println("\n\n", "#"^84)
    @printf("PASS 3 — every node covers X: ΓX = %.4f, %d slices\n", Γ3, N3)
    println("#"^84)
    for nd in nodes
        @printf(
            "  piece %-6s covers X at level %.4f\n",
            string(nd),
            PCQ.gamma_cover_set(pclf.pieces[nd], Xh)
        )
    end
    q1 = arm("ROUTE 1 — every-node level", common; ΓX = Γ3, nb_levels = N3)
    q2 = arm("ROUTE 2 — every-node level", pclf; ΓX = Γ3, nb_levels = N3)
end

function side_by_side(title, a, b)
    println("\n", "="^84)
    println(title)
    println("="^84)
    rows = (
        ("build time (s)", med(a.t), med(b.t), true),
        ("cells", a.s.cells, b.s.cells, true),
        ("transitions", a.s.transitions, b.s.transitions, true),
        ("slices", a.s.slices, b.s.slices, false),
        ("nodes carrying a partition", a.s.nodes, b.s.nodes, false),
        ("Σ parts over all cells", a.s.sum_p, b.s.sum_p, true),
        ("mean parts / cell", a.s.mean_p, b.s.mean_p, true),
        ("median parts / cell", a.s.med_p, b.s.med_p, false),
        ("MAX parts / cell", a.s.max_p, b.s.max_p, true),
        ("Σ facets over all cells", a.s.sum_f, b.s.sum_f, true),
        ("mean facets / cell", a.s.mean_f, b.s.mean_f, true),
        ("median facets / cell", a.s.med_f, b.s.med_f, false),
        ("MAX facets / cell", a.s.max_f, b.s.max_f, true),
    )
    @printf("%-28s %14s %14s %12s\n", "", "route 1", "route 2", "1 / 2")
    println("-"^84)
    for (nm, x, y, ratio) in rows
        @printf(
            "%-28s %14.2f %14.2f %12s\n",
            nm,
            x,
            y,
            ratio ? @sprintf("%.2f", x / y) : "—"
        )
    end
    println("="^84)
    @printf(
        "τX: route 1 %.4f, route 2 %.4f — %s\n",
        maximum(a.τ),
        maximum(b.τ),
        isapprox(maximum(a.τ), maximum(b.τ); rtol = 0.02) ?
        "within 2 %, the arms tile the same region and the ratios above stand" :
        "MORE THAN 2 % APART: the arms tile different regions, so these ratios measure COVERAGE",
    )
    return nothing
end

side_by_side("PASS 1 — the poster's settings (each arm its own ladder)", r1, r2)
p1 === nothing || side_by_side("PASS 2 — same ladder, the level the UNION needs", p1, p2)
q1 === nothing || side_by_side("PASS 3 — same ladder, the level EVERY NODE needs", q1, q2)

println(
    """

HOW TO READ THIS BEFORE ADOPTING THE EXAMPLE

  PASS 1 is the poster's configuration and it is NOT a cost comparison: the arms derive different
  outer levels from different gauges, so route 2 tiles a larger region with one more slice and its
  cell count reflects that, not its efficiency.

  PASS 2 is the cost comparison. Same ΓX, same slice count, same tolerance, same certificate: the
  only difference left is whether the PCLF was determinised first.

  The per-cell rows -- mean and MAX parts and facets -- are meaningful in BOTH passes, because they
  do not depend on how many cells an arm built. They are the mechanism: route 1's cells inherit the
  induced common's fragmentation and route 2's do not.""",
)
