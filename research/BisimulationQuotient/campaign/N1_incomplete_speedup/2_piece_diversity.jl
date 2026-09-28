# Finding an instance whose PCLF pieces are genuinely DIFFERENT, and checking the benefit survives.
#
# Why this matters. In the measured instance the two node pieces are the same square to within a
# hair, and that is not something the certificate discovered -- it is imposed:
#
#     rotation_templates(nodes; mode = :rotation)  ->  Dict(v => Rot for v in nodes)
#
# every node gets the SAME normal directions `G`, so the pieces {x : -w <= G x <= w} can differ only
# in the offsets `w`. With four facets and a shared orientation the two nodes have no shape freedom
# left. A reader is then entitled to object that the PCLF is a common Lyapunov function wearing two
# hats, and that route 2's advantage is an artefact of that degeneracy rather than a property of the
# construction. The objection is reasonable and this script exists to answer it.
#
# Three independent levers, swept here:
#
#   ORIENTATION    `mode = :alternating` gives node (1,) the identity template and node (2,) `Rot`,
#                  so the two pieces are rotated by theta with respect to one another BY CONSTRUCTION.
#   TRANSVERSALITY the shear family `A_i = rho * [1 ±alpha; 0 1] * D`. At alpha = 0 the modes coincide
#                  and nothing can force the pieces apart; as alpha grows the dynamics demand it.
#   COMPACTNESS   more facets, via `compute_polyhedral_pieces_pclf` with a shared conic partition, so
#                  the solver can DISCOVER a node-dependent shape rather than have one imposed. Swept
#                  in PART TWO below, where it is used to ask whether path-completeness buys a RATE.
#
# What is reported, and why these columns. `piece gap` is the largest relative difference between the
# two pieces' support functions over 64 directions: 0 means identical, and it is exact for convex sets
# and needs no volume computation. Its CONSEQUENCE is the shape of the induced common, which is what
# route 1 actually pays for -- so `parts` and `facets` of that common are reported beside it. On the
# dual graph the common is the intersection, so 4 facets means one piece contains the other and 8 means
# they genuinely cross; on the primal it is the union, so 1 part means nesting and 3 or more means the
# pieces cross.
#
# THE HONEST READING, STATED IN ADVANCE. Pieces that genuinely differ will make route 1's induced
# common more complex, and therefore make route 2 look BETTER. That is not a favourable tuning choice:
# it is the removal of a confound, since the shared template was suppressing a difference the dynamics
# would otherwise produce. Both the degenerate and the diverse instance are reported so the effect of
# removing it is visible rather than absorbed into a headline number.
#
# ------------------------------------------------------------------------------------------------
# TWO THINGS THE FIRST RUN OF THIS SCRIPT GOT WRONG, KEPT HERE SO THEY ARE NOT REPEATED
# ------------------------------------------------------------------------------------------------
#
# 1. THE LADDER WAS NOT PINNED, so the two arms were compared at different depths -- 6 rungs against
#    7, 8 against 9 -- and one rung is worth 4-5x in cells. The ratios it produced were one rung of
#    growth and nothing else, and the single configuration where route 2 "won" was the single one
#    where both arms happened to land on the same depth. Region-free plus `ΓX` and `nb_levels` is now
#    enforced, and an unequal rung count raises a warning rather than passing silently.
#
# 2. `alpha = 0` REPORTED "no certificate", WRONGLY -- SINCE FIXED IN THE LIBRARY. At alpha = 0 both modes equal
#    A = diag(0.6, 0.48), which is trivially stable. With the orthogonal template G = R(pi/6) the
#    condition is governed by  G A Gᵀ = [0.57 0.052; 0.052 0.51],  whose induced infinity-norm is
#    0.622 -- a perfectly good certificate. But the bisection brackets from above by
#    `b = max_i ||A_i||_inf = 0.6 < 0.622`, so the TOP of the bracket is itself infeasible, no trial
#    rate ever passes, and `Inf` comes back. At alpha = 0.2 the bracket is 0.696 > 0.622 and the same
#    code succeeds, which is the whole of the discontinuity.
#
#    This is a limitation of the bracket, not a fact about the system: `max_i ||A_i||_inf` is only a
#    valid upper bound when the template is aligned with the norm that bound is taken in, and a
#    rotated template is not. It is left unfixed here (alpha = 0 is the degenerate endpoint where the
#    two modes coincide and the dual graph can have nothing to gain), but it will mis-report any
#    rotated-template problem whose true rate exceeds the norm bound, and it belongs in the library's
#    known issues rather than in one script's comments.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Printf
using LinearAlgebra

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)

"""
Two modes sheared in OPPOSITE directions by `α`: `A_i = ρ · [1 ±α; 0 1] · D`.

Shear, not rotation, and the reason is a trap worth recording. The first version used
`A_i = ρ R(±α) D`, which produced **no certificate at any α**: the template is a fixed-orientation
4-facet square, a rotation turns it, and fitting a turned square inside a scaled copy of itself costs
a gauge factor up to √2, so `γ ≥ ρ√2 = 1.06 > 1`. Rotations are exactly the case a single fixed
template cannot follow, and that is a *feasibility* separation — a different experiment, since a cost
comparison needs both routes to run.

A shear leaves the square a parallelogram, which fits in a scaled square at roughly `1 + |α|`, so
`γ ≈ ρ(1 + α)` stays under 1 across the sweep. The two modes remain distinct and are never negatives
of each other (`A₊ + A₋ = 2ρD ≠ 0`), so the symmetric-template degeneracy cannot reappear either.
"""
function transverse_pair(α; ρ = 0.6, D = [1.0 0.0; 0.0 0.8])
    shear(s) = [1.0 s; 0.0 1.0]
    return [ρ * shear(α) * D, ρ * shear(-α) * D]
end

"""
Number of rungs both arms are forced onto. Matches the main comparison's region-free case.

**Without this the screen measures nothing.** A first version left the ladder free and the two arms
came out at *different depths* — 6 rungs against 7 at α = 0.4, 8 against 9 at α = 0.6 — and a single
extra rung multiplies the partition by 4–5. The cell ratios it produced (0.25, 0.32) were almost
exactly one rung of growth, and the one configuration where route 2 "won" was the one configuration
where the two arms happened to land on the same depth. The measurement was of the ladder, not of the
routes. `plan.md` §2d says exactly this, and the first version of this script ignored it.
"""
const NB_LEVELS = 5

"""
Region-FREE, so the only things creating cells are the graph and the matrices.

Observation regions are dropped here on purpose. They would add a third cause of cell creation on top
of the two this screen is about, and they are also what makes the ladder depth arm-dependent (it
descends until it clears them). Region-free plus a pinned ladder is the configuration `plan.md` §2d
established as the comparable one.
"""
function problem_at(α)
    A = transverse_pair(α)
    f = ST.with_switching(
        HybridSystems.discreteswitchedsystem(A),
        HybridSystems.ControlledSwitching(),
    )
    X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
    return (; f, A, problem = PR.BisimulationQuotientProblem(f, X, typeof(X)[]))
end

function oracle_pclf(f, dual, tmode; θ = π / 6)
    graph = PCLF.generate_DeBruijn_edges(2, 1; dual = dual)
    nodes = sort(collect(graph.verts); by = string)
    return PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
        f,
        graph,
        OPTIMIZER;
        Gmats = rotation_templates(nodes; θ = θ, mode = tmode),
        MLF = true,
        verbose = false,
    )
end

"""
How different two convex pieces are: the largest relative gap between their support functions.

Exact for convex sets, needs no volume, and reads as a percentage — `0.00` means the two pieces are
the same set, `0.20` that some direction sees one extend 20 % further than the other. Volume would
answer a different and weaker question: two sets can have equal area and different shape.
"""
function piece_gap(P1, P2; K::Int = 64)
    worst = 0.0
    for k in 0:(K - 1)
        d = [cos(2π * k / K), sin(2π * k / K)]
        h1, h2 = LazySets.ρ(d, P1), LazySets.ρ(d, P2)
        worst = max(worst, abs(h1 - h2) / max(abs(h1), abs(h2), eps()))
    end
    return worst
end

"The shape of the induced common at γ = 1: how many convex parts, and the widest of them."
function common_shape(common)
    S = PCLF.get_sublevel_set(common.pieces[:clf], 1.0; atol = 1e-6)
    parts = S isa LazySets.UnionSetArray ? S.array : [S]
    return (;
        nparts = length(parts),
        nfacets = maximum(length(LazySets.constraints_list(P)) for P in parts),
    )
end

"One quotient on the pinned ladder: same rung count and same outer level for both arms."
function pinned(problem, pclf, ΓX)
    q = build_quotient(
        problem,
        pclf;
        atol = 1e-3,
        nb_levels = NB_LEVELS,
        ΓX = ΓX,
        max_slices = NB_LEVELS,
        print_level = 0,
    ).quotient
    return (; cells = length(q.states), slices = PCQ.num_slices(q))
end

"""
Median build time of two arms, sampled ALTERNATELY, with route 1's determinisation inside its own
closure.

Cell count settles nothing on a complete graph: lifting is redundant there, so route 2 pays `|S|`
copies and the ratio sits near `1/|S|` whatever else is true. What the primal graph decides is
**time**, because route 1's cells fragment into the induced union's parts. A screen that reports only
counts therefore cannot answer "does the benefit survive diverse pieces", which is the whole question.
"""
function timed_pair(problem, p1, p2, ΓX, make_p1; samples = 3)
    b1 = () -> (make_p1(); pinned(problem, p1, ΓX))
    b2 = () -> pinned(problem, p2, ΓX)
    b1()
    b2()
    t1, t2 = Float64[], Float64[]
    for _ in 1:samples, (b, ts) in ((b1, t1), (b2, t2))
        GC.gc()
        t0 = time_ns()
        b()
        push!(ts, (time_ns() - t0) / 1e9)
    end
    med(v) = sort(v)[(length(v) + 1) ÷ 2]
    return (;
        m1 = med(t1),
        lo1 = minimum(t1),
        hi1 = maximum(t1),
        m2 = med(t2),
        lo2 = minimum(t2),
        hi2 = maximum(t2),
    )
end

function measure(α, dual, tmode)
    (; f, problem) = problem_at(α)
    p2 = oracle_pclf(f, dual, tmode)
    (isfinite(p2.JSRapprox) && p2.JSRapprox < 1) ||
        return (; ok = false, reason = @sprintf("no certificate (rate %.4f)", p2.JSRapprox))
    p1 = PCLF.build_common_lyapunov(p2)
    isapprox(p1.JSRapprox, p2.JSRapprox; rtol = 1e-4) ||
        return (; ok = false, reason = "rates differ — not a cost comparison")

    nodes = sort(collect(p2.graph.verts); by = string)
    pieces = [PCLF.get_sublevel_set(p2.pieces[v], 1.0) for v in nodes]
    shape = common_shape(p1)

    # ΓX comes from route 1's own natural outer level, so the baseline keeps the geometry it would
    # have chosen for itself and it is route 2 that is made to conform.
    probe = build_quotient(
        problem,
        p1;
        atol = 1e-3,
        nb_levels = NB_LEVELS,
        max_slices = NB_LEVELS,
        print_level = 0,
    )
    ΓX = maximum(MOI.get(probe.optimizer, MOI.RawOptimizerAttribute("Γ")))
    r1, r2 = pinned(problem, p1, ΓX), pinned(problem, p2, ΓX)
    return (;
        ok = true,
        rate = p2.JSRapprox,
        gap = piece_gap(pieces[1], pieces[2]),
        shape.nparts,
        shape.nfacets,
        c1 = r1.cells,
        c2 = r2.cells,
        s1 = r1.slices,
        s2 = r2.slices,
        ratio = r1.cells / r2.cells,
    )
end

println(
    """
Do the pieces have to be near-identical, and does the benefit survive making them differ?

`piece gap` is the largest relative difference between the two node pieces' support functions: 0.00
means the shared template has forced them to coincide. `parts`/`facets` describe the induced common
that route 1 must then work on -- on the DUAL graph that is the intersection (4 facets = one piece
inside the other, 8 = they genuinely cross), on the PRIMAL the union (1 part = nesting, 3+ = crossing).

`cells r1/r2` above 1 means route 2 builds fewer.
""",
)

@printf(
    "%-12s %-7s %-6s %-8s %-10s %-7s %-8s %-9s %-9s %-8s %s\n",
    "template",
    "graph",
    "alpha",
    "rate",
    "piece gap",
    "parts",
    "facets",
    "r1 cells",
    "r2 cells",
    "rungs",
    "cells r1/r2"
)
println("-"^112)

for tmode in (:rotation, :alternating), dual in (true, false)
    for α in (0.0, 0.2, 0.4, 0.6, 0.8)
        r = measure(α, dual, tmode)
        if !r.ok
            @printf(
                "%-12s %-7s %-6.2f %s\n",
                string(tmode),
                dual ? "dual" : "primal",
                α,
                r.reason
            )
            continue
        end
        # The pinning has to be CHECKED, not assumed: an unequal rung count is what invalidated the
        # first version of this screen, and it is silent unless something looks.
        r.s1 == r.s2 == NB_LEVELS || @warn(
            "ladder not pinned — this row compares different depths and must be discarded",
            tmode,
            graph = dual ? "dual" : "primal",
            α,
            rungs_route1 = r.s1,
            rungs_route2 = r.s2
        )
        @printf(
            "%-12s %-7s %-6.2f %-8.5f %-10.3f %-7d %-8d %-9d %-9d %-8s %.2f\n",
            string(tmode),
            dual ? "dual" : "primal",
            α,
            r.rate,
            r.gap,
            r.nparts,
            r.nfacets,
            r.c1,
            r.c2,
            r.s1 == r.s2 ? string(r.s1) : "$(r.s1)/$(r.s2) !",
            r.ratio
        )
    end
    println()
end

# -------------------------------------------------------------------------------------------------
# PART TWO: the FUNCTIONAL gain — a rate no common Lyapunov function reaches
# -------------------------------------------------------------------------------------------------
#
# Part one asked whether the pieces CAN differ. This part asks the sharper question: is there an
# instance where path-completeness does work no common Lyapunov function could do, so that the
# node-dependent pieces are a necessity rather than a stylistic choice?
#
# ONE STRUCTURAL FACT DECIDES WHERE TO LOOK, AND IT RULES OUT THE DUAL GRAPH ENTIRELY.
#
# On a co-complete graph — the dual De Bruijn is one — the induced common is `V_max = max_i V_i`, a
# maximum of convex gauges, hence ITSELF CONVEX, and it certifies the PCLF's own rate. So whatever a
# dual-graph PCLF certifies, a convex common Lyapunov function certifies too: the dual graph can never
# give a functional advantage over CLFs, only a computational one. That is not a flaw — it is what
# makes the dual graph the clean instance for isolating the ALGORITHMIC gain, because no functional
# gain can be confounded with it.
#
# On a complete graph — the primal De Bruijn — the induced common is `V_min = min_i V_i`, a minimum of
# convex gauges and therefore NON-CONVEX. Only there can a PCLF certify a rate outside the reach of
# convex common functions, and only there must the pieces genuinely differ.
#
# THE TEST. At each conic-partition order, compare
#
#     the SINGLE NODE (k = 0)   — literally a common Lyapunov function in that template class
#     the PRIMAL De Bruijn      — path-complete, node-dependent
#
# A functional gain is `clf rate >= 1 > pclf rate`: the common function cannot certify stability and
# the path-complete one can. Where that happens the pieces cannot be interchangeable, and the measured
# piece gap should say so.
#
# An anisotropic family was tried first and rejected BEFORE running it: with
# `A_i = ρ R(±φ) diag(1, s) R(∓φ)` every `A_i` is symmetric, so `‖A_i‖₂ = ρ < 1` and the Euclidean ball
# is already a common Lyapunov function. Such a family can exhibit piece diversity but never a
# functional gain, because there is nothing left for path-completeness to buy.

"""
The two-mode pair from the observer-graph study, optionally scaled.

Both modes are individually stable (spectral radii 0.80 and 0.45) but their products are not tame, so
this is the standard shape of a system where a single common function struggles and a path-complete
one need not. `scale` walks it towards the edge of certifiability, which is where the gap between the
two is widest.
"""
function hard_pair(; scale = 1.0)
    A1 = scale * (1.0 / 10.0) * [1.5519 0.4474; 7.6412 7.4716]
    A2 = scale * (1.0 / 10.0) * [0.4750 9.1755; 1.8955 0.1850]
    return [A1, A2]
end

function hard_problem(scale)
    f = ST.with_switching(
        HybridSystems.discreteswitchedsystem(hard_pair(; scale = scale)),
        HybridSystems.ControlledSwitching(),
    )
    X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
    return (; f, problem = PR.BisimulationQuotientProblem(f, X, typeof(X)[]))
end

"A PCLF with `order`-refined conic pieces, the partition SHARED across nodes so nothing is imposed."
function conic_pclf(f, graph, order)
    nodes = sort(collect(graph.verts); by = string)
    return PCLF.compute_polyhedral_pieces_pclf(
        f,
        graph,
        OPTIMIZER,
        PCLF.conic_partitions_dict_2d(order, nodes);
        MLF = true,
    )
end

function measure_functional(scale, order)
    (; f, problem) = hard_problem(scale)

    # k = 0 is one node carrying a self-loop per mode: a COMMON Lyapunov function in this class.
    clf = conic_pclf(f, PCLF.generate_DeBruijn_edges(2, 0), order)
    p2 = conic_pclf(f, PCLF.generate_DeBruijn_edges(2, 1), order)
    row = (; scale, order, clf_rate = clf.JSRapprox, pclf_rate = p2.JSRapprox)

    # A PCLF can never be WORSE than the common function: setting every node's piece to the common
    # one is a feasible point of the richer problem. So `pclf > clf` is not a fact about the system,
    # it is the solver failing to converge on the bigger program, and reporting it as "the PCLF does
    # not certify" would be reporting a numerical failure as a result.
    p2.JSRapprox <= clf.JSRapprox * (1 + 1e-6) || return (;
        row...,
        ok = false,
        reason = "SOLVER FAILURE: pclf rate exceeds clf rate, which is impossible — discard",
    )

    (isfinite(p2.JSRapprox) && p2.JSRapprox < 1) ||
        return (; row..., ok = false, reason = "the PCLF does not certify either")
    p1 = PCLF.build_common_lyapunov(p2)
    isapprox(p1.JSRapprox, p2.JSRapprox; rtol = 1e-4) ||
        return (; row..., ok = false, reason = "rates differ — not a cost comparison")

    nodes = sort(collect(p2.graph.verts); by = string)
    pieces = [PCLF.get_sublevel_set(p2.pieces[v], 1.0) for v in nodes]
    shape = common_shape(p1)
    probe = build_quotient(
        problem,
        p1;
        atol = 1e-3,
        nb_levels = NB_LEVELS,
        max_slices = NB_LEVELS,
        print_level = 0,
    )
    ΓX = maximum(MOI.get(probe.optimizer, MOI.RawOptimizerAttribute("Γ")))
    r1, r2 = pinned(problem, p1, ΓX), pinned(problem, p2, ΓX)
    t = timed_pair(problem, p1, p2, ΓX, () -> PCLF.build_common_lyapunov(p2))
    return (;
        row...,
        ok = true,
        gap = piece_gap(pieces[1], pieces[2]),
        nparts = shape.nparts,
        c1 = r1.cells,
        c2 = r2.cells,
        s1 = r1.slices,
        s2 = r2.slices,
        ratio = r1.cells / r2.cells,
        t,
        tratio = t.m1 / t.m2,
        disjoint = t.lo1 > t.hi2 || t.lo2 > t.hi1,
    )
end

println(
    """

================================================================================================
PART TWO — the functional gain: a rate no common Lyapunov function reaches
================================================================================================
`clf rate` is the SINGLE-NODE graph, which IS a common Lyapunov function in this template class.
`pclf rate` is the primal De Bruijn of order 1. A functional gain is `clf >= 1 > pclf`: the common
function cannot certify stability and the path-complete one can, so the pieces cannot be
interchangeable.

The dual graph is deliberately absent. Its induced common is max_i V_i, convex, certifying the same
rate — so no functional gain is possible there by construction, and that is precisely what makes it
the clean instance for the purely algorithmic claim.
""",
)

@printf(
    "%-7s %-7s %-10s %-10s %-9s %-10s %-7s %-9s %-11s %s\n",
    "scale",
    "order",
    "clf rate",
    "pclf rate",
    "gain?",
    "piece gap",
    "parts",
    "cells r1/r2",
    "TIME r1/r2",
    "verdict?"
)
println("-"^108)

for scale in (1.0, 1.05, 1.1), order in (1, 2, 3)
    r = measure_functional(scale, order)
    gain = isfinite(r.pclf_rate) && r.pclf_rate < 1 && !(r.clf_rate < 1)
    if !r.ok
        @printf(
            "%-7.2f %-7d %-10.4f %-10.4f %-9s %s\n",
            r.scale,
            r.order,
            r.clf_rate,
            r.pclf_rate,
            gain ? "YES" : "no",
            r.reason
        )
        continue
    end
    r.s1 == r.s2 == NB_LEVELS ||
        @warn("ladder not pinned — discard this row", r.scale, r.order, r.s1, r.s2)
    @printf(
        "%-7.2f %-7d %-10.4f %-10.4f %-9s %-10.3f %-7d %-11.2f %-11.2f %s\n",
        r.scale,
        r.order,
        r.clf_rate,
        r.pclf_rate,
        gain ? "*** YES" : "no",
        r.gap,
        r.nparts,
        r.ratio,
        r.tratio,
        r.disjoint ? "disjoint" : "OVERLAP — tie"
    )
end

println(
    """

Reading it.

A row marked `*** YES` is the paper's second example: the common Lyapunov function in that template
class does NOT certify the system and the path-complete one does, so the node-dependent pieces carry
real information instead of repeating one shape. `piece gap` should then be far from zero and `parts`
above 1 — the induced common is a genuinely non-convex union, which is exactly what a convex common
function cannot be.

`cells r1/r2` below 1 is EXPECTED on the primal graph and is not a failure: the graph is complete, so
lifting is redundant and route 2 pays |S| copies (plan.md §1). The primal graph's payoff is time, not
count — route 1's cells fragment into the union's parts (plan.md §2c). A functional example and a cell
advantage are different claims, and this table settles only the first.

If no row qualifies, the honest conclusion is that within reach of these templates path-completeness
buys no rate a common function cannot match, and the paper's claim should stay computational.""",
)
