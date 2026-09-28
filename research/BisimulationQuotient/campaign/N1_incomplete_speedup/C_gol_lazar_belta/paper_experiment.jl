# GLB, primal De Bruijn — the experiment intended for the paper.
#
#     julia --project=test .../C_gol_lazar_belta/paper_experiment.jl        SAMPLES=n  FIGURES=0
#
# Everything is theirs except the certificate: their dynamics, their three observation regions, their
# co-safe formula, their initial point. The certificate is a path-complete Lyapunov function from an
# oracle, and the two routes are compared on it.
#
# WHAT THE DESIGN FIXES, AND WHY EACH CHOICE IS FORCED RATHER THAN TUNED
#
#   solver          HiGHS, not Clarabel. These constraints are all linear, so this is a pure LP.
#                   `compute_polyhedral_pieces_pclf` treats any termination status other than
#                   OPTIMAL/FEASIBLE_POINT as infeasible, and Clarabel — an interior-point conic
#                   solver — stops returning OPTIMAL as the LP grows and flattens. The effect is a
#                   rate that gets WORSE as the template gets richer (0.8674 at conic order 2,
#                   0.9575 at order 3, 0.9671 at order 4), which is impossible: refining a conic
#                   partition can only enlarge the feasible set, since duplicating each row onto the
#                   two half-cones reproduces the coarser solution. With HiGHS the sequence is
#                   monotone (0.9381, 0.8674, 0.8596, 0.8577) and converges to 0.8558, the largest
#                   spectral radius over the self-loops — the JSR lower bound. That is the expected
#                   behaviour and the reason to trust it.
#
#   template        SEVEN PLACED CONES, not a refinement order. `conic_partitions_2d(k)` doubles the
#                   cone count at every level — 4, 8, 16, 32 — and the template row count is the cone
#                   count, so every cell of both arms pays for the refinement. Placing the rays
#                   instead, by coordinate descent on the certified rate, dominates that:
#
#                     rows   rate       parts   slices
#                        8   0.867381       5        7      repo order 2
#                       16   0.859586      11        7      repo order 3
#                        7   0.864529      11        6      <- used here
#                        8   0.859099      11        6
#
#                   Seven placed rays beat eight uniform ones AND the repo's order 2, need one slice
#                   fewer than either repo order, and fragment the induced common as much as order 3
#                   does with half its rows. Fewer rows, fewer slices, same fragmentation.
#
#   scale           the pieces are rescaled by a COMMON factor so the region-containment level is 6.
#                   A global scale leaves the edge conditions and the rate untouched, but the
#                   clearing test downstream uses an ABSOLUTE tolerance, and HiGHS returns a solution
#                   some 500x smaller than Clarabel's — at that scale the test silently answers
#                   nonsense.
#
#   ΓX              Γ*, the smallest level at which EVERY node's piece contains EVERY region. Not
#                   the level at which their union does: the universal semantics quantifies over all
#                   nodes, and a node whose piece misses R3 cannot certify anything there. Pinning at
#                   the union level left exactly that hole in node (2).
#
#   slices          the construction's own stopping rule — the first level whose terminal set no
#                   longer meets R1, R2 or R3. Forcing fewer leaves D overlapping a region, which
#                   mislabels every point of the overlap.
#
# ONE THING THIS EXPERIMENT CANNOT CLAIM. The LP has many optima and the solver picks one. At conic
# order 2 and the same rate to five decimals, Clarabel's optimum induces a common of 17 parts and
# HiGHS's one of 5. Fragmentation is the mechanism of route 2's advantage, so its MAGNITUDE is a
# property of the certificate returned, not of the problem. The claim that survives is qualitative:
# determinising a complete graph yields a non-convex union and route 1's cells inherit it. HiGHS is
# the conservative choice here — it fragments less, so it understates our margin.

include(joinpath(dirname(dirname(dirname(@__DIR__))), "common.jl"))

using Printf
using Spot
import Statistics
import HiGHS

gr()

const SAMPLES = parse(Int, get(ENV, "SAMPLES", "1"))
const FIGURES = get(ENV, "FIGURES", "1") == "1"
const ATOL = 1e-3          # the campaign's tolerance; 1e-4 is finer and several times slower
const SCALE = 6.0          # the level the region-containment requirement is normalised onto

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
const Rs = (("R1", R1), ("R2", R2), ("R3", R3))
const REGIONS = [UT._as_hpolytope(R) for (_, R) in Rs]
φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)   # point `a` of the paper

graph = PCLF.generate_DeBruijn_edges(2, 1)      # primal De Bruijn, order 1
nodes = sort(collect(graph.verts); by = string)

"The smallest level at which every node's piece contains every observation region."
Γstar(p) = maximum(
    PCQ.gamma_cover_set(p.pieces[nd], UT._as_hpolytope(R)) for nd in nodes for (_, R) in Rs
)

function normalise!(p, target)
    α = Γstar(p) / target
    for nd in keys(p.pieces)
        p.pieces[nd].G = p.pieces[nd].G / α
    end
    return p
end

# The ray angles below are not a refinement level — they were found by coordinate descent on the
# certified rate, starting from the uniform 7-ray partition. See the header for why that matters.
const RAYS_DEG = [0.0, 32.6, 65.8, 91.5, 112.6, 131.4, 152.0, 180.0]
const CONES = let a = deg2rad.(RAYS_DEG)
    [hcat([cos(a[i]), sin(a[i])], [cos(a[i + 1]), sin(a[i + 1])]) for i in 1:(length(a) - 1)]
end

println("computing the certificate (HiGHS, 7 placed cones) ...")
pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    graph,
    HiGHS.Optimizer,
    Dict(nd => CONES for nd in nodes);
    MLF = true,
)
isinf(pclf.JSRapprox) && error("no certificate found")
normalise!(pclf, SCALE)
common = PCLF.build_common_lyapunov(pclf)
γ = pclf.JSRapprox
ΓX = Γstar(pclf)

NB = something(
    findfirst(
        j ->
            PCQ.all_nodes_clear_regions(pclf, ΓX * γ^(j - 1), REGIONS; tol = 1e-2) &&
            PCQ.all_nodes_clear_regions(common, ΓX * γ^(j - 1), REGIONS; tol = 1e-2),
        1:80,
    ),
    -1,
)
NB > 0 || error("no region-free terminal level within 80 rungs")

common_parts = let S = PCLF.get_sublevel_set(common.pieces[:clf], 1.0; atol = 1e-6)
    S isa LazySets.UnionSetArray ? length(S.array) : 1
end

println("\n", "="^84)
println("THE CERTIFICATE")
println("="^84)
@printf(
    "graph            primal De Bruijn order 1 — %d nodes, %d edges, complete: %s\n",
    length(graph.verts),
    length(graph.edges),
    PCLF.is_complete(graph, 1:2)
)
@printf(
    "template         %d placed cones — %d rows per piece, rays at %s deg\n",
    length(CONES),
    size(pclf.pieces[nodes[1]].G, 1),
    join([@sprintf("%.1f", d) for d in RAYS_DEG], " ")
)
@printf("certified rate   %.6f\n", γ)
@printf("induced common   %d convex parts — what route 1 builds on\n", common_parts)
@printf("ΓX = Γ*          %.4f   (every node contains every region)\n", ΓX)
@printf("slices           %d, τD = %.4f\n", NB, ΓX * γ^(NB - 1))
println("\nregion containment per node, at ΓX:")
for nd in nodes, (nm, R) in Rs
    S = PCLF.get_sublevel_set(pclf.pieces[nd], ΓX)
    ok = try
        UT.is_included(UT._as_hpolytope(R), S)
    catch
        missing
    end
    @printf(
        "  node %-6s %-3s : %s\n",
        string(nd),
        nm,
        ok === true ? "contained" : "NOT CONTAINED"
    )
end

build(cert) = build_quotient(
    problem,
    cert;
    atol = ATOL,
    nb_levels = NB,
    ΓX = ΓX,
    max_slices = NB,
    print_level = 0,
)

function stats(q)
    parts, faces = PCQ.cell_complexities(q)
    st = PCQ.bisimulation_stats(q)
    return (;
        cells = st[:num_states],
        transitions = st[:num_transitions],
        sum_p = sum(parts),
        mean_p = Statistics.mean(parts),
        max_p = maximum(parts),
        sum_f = sum(faces),
        mean_f = Statistics.mean(faces),
        max_f = maximum(faces),
    )
end

function arm(name, cert)
    println("\n", "-"^84, "\n", name, "\n", "-"^84)
    print("  warming up ... ")
    b = build(cert)
    println("done")
    ts = Float64[]
    for k in 1:SAMPLES
        GC.gc()
        t0 = time_ns()
        b = build(cert)
        push!(ts, (time_ns() - t0) / 1e9)
        @printf("  round %d: %.3f s\n", k, ts[k])
    end
    Dp = b.D isa LazySets.UnionSetArray ? b.D.array : [UT._as_hpolytope(b.D)]
    all(LazySets.isdisjoint(P, R) for P in Dp, R in REGIONS) ||
        error("[$name] terminal set D meets an observation region")
    @printf("  D: %d part(s), clears every region\n", length(Dp))
    return (; t = ts, s = stats(b.quotient), q = b.quotient, D = b.D)
end

r1 = arm("ROUTE 1 — determinise first, then build on the induced common", common)
r2 = arm("ROUTE 2 — build directly on the PCLF, one partition per node", pclf)

# ── the specification, both arms, both quantifiers ───────────────────────────────────────────────
f_forall = ST.with_switching(f, HybridSystems.AutonomousSwitching())
function solve_timed(f_sys, a, label)
    run() = synthesize_cosafe_ltl(
        f_sys,
        a.q,
        Dionysos.spot_stepper(φ),
        Dict{Symbol, Any}(:D => a.D, :R1 => R1, :R2 => R2, :R3 => R3),
        Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
        x0;
        print_level = 1,
    )
    run()
    ts = Float64[]
    for _ in 1:SAMPLES
        GC.gc()
        t0 = time_ns()
        run()
        push!(ts, (time_ns() - t0) / 1e9)
    end
    res = run()
    @printf(
        "  %-22s %8.3f s   %d cells\n",
        label,
        Statistics.median(ts),
        length(res.controllable_set)
    )
    return (; res, ts)
end

println("\n", "-"^84, "\nTHE SPECIFICATION\n", "-"^84)
sols = Dict(
    "route 1 ∃" => solve_timed(f, r1, "route 1 synthesis ∃"),
    "route 1 ∀" => solve_timed(f_forall, r1, "route 1 verification ∀"),
    "route 2 ∃" => solve_timed(f, r2, "route 2 synthesis ∃"),
    "route 2 ∀" => solve_timed(f_forall, r2, "route 2 verification ∀"),
)

med(v) = Statistics.median(v)
println("\n", "="^84)
println("RESULT")
println("="^84)
@printf("%-26s %14s %14s %10s\n", "", "route 1", "route 2", "1 / 2")
println("-"^84)
for (nm, a, b) in (
    ("build (s)", med(r1.t), med(r2.t)),
    ("synthesis ∃ (s)", med(sols["route 1 ∃"].ts), med(sols["route 2 ∃"].ts)),
    ("verification ∀ (s)", med(sols["route 1 ∀"].ts), med(sols["route 2 ∀"].ts)),
    (
        "TOTAL (s)",
        med(r1.t) + med(sols["route 1 ∃"].ts) + med(sols["route 1 ∀"].ts),
        med(r2.t) + med(sols["route 2 ∃"].ts) + med(sols["route 2 ∀"].ts),
    ),
    ("cells", r1.s.cells, r2.s.cells),
    ("Σ parts", r1.s.sum_p, r2.s.sum_p),
    ("mean parts / cell", r1.s.mean_p, r2.s.mean_p),
    ("MAX parts / cell", r1.s.max_p, r2.s.max_p),
    ("mean facets / cell", r1.s.mean_f, r2.s.mean_f),
    ("MAX facets / cell", r1.s.max_f, r2.s.max_f),
)
    @printf("%-26s %14.2f %14.2f %10.2f\n", nm, a, b, a / b)
end
println("="^84)

if FIGURES
    save(fig, n) = (savefig(fig, joinpath(@__DIR__, n * ".png")); println("wrote $n.png"))
    save(
        plot_synthesis_result(
            r1.q,
            sols["route 1 ∃"].res,
            problem,
            "Route 1 — synthesis ∃",
        ),
        "paper_route1_synthesis",
    )
    save(
        plot_synthesis_result(
            r1.q,
            sols["route 1 ∀"].res,
            problem,
            "Route 1 — verification ∀",
        ),
        "paper_route1_verification",
    )
    for (tag, key) in (("synthesis", "route 2 ∃"), ("verification", "route 2 ∀"))
        res = sols[key].res
        panels = map(quotient_nodes(r2.q)) do nd
            ids(s) = [i for i in s if r2.q.states[i].node == nd]
            p = plot(;
                aspect_ratio = :equal,
                legend = false,
                title = "Route 2 — $tag, node $nd",
                titlefontsize = 9,
            )
            _plot_winning_losing!(
                p,
                r2.q,
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
        save(
            plot(
                panels...;
                layout = (1, length(panels)),
                size = panel_row_size(lims, length(panels)),
            ),
            "paper_route2_$tag",
        )
    end
end
