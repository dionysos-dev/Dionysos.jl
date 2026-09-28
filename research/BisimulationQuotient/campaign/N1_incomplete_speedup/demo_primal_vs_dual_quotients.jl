# THE FOUR QUOTIENTS IN ONE PICTURE — {route 1, route 2} × {primal, dual}, same problem throughout.
#
#     julia --project=test .../demo_primal_vs_dual_quotients.jl [glb|A] [rungs]      SAMPLES=n
#
# The claim it is built to test — NOT one it is built to confirm: route 2 beats route 1 on both
# orientations of the graph the oracle hands us, for two DIFFERENT reasons. It holds on the canonical
# example and it FAILED on GLB's dual row in every measurement so far, which is why the conditional
# form is the one to quote: route 2 wins on a complete graph via cell simplicity, and on a
# co-complete graph only when the per-node reduction `g` exceeds |S|.
#
# Each row is timed over `SAMPLES` paired rounds (default 3) with both arms run back-to-back inside a
# round, after a warm-up call that absorbs compilation, and every figure is drawn strictly after the
# last measurement — drawing between timed builds was a real defect in an earlier version of this
# campaign.
#
# WHAT THE GRID SHOWS
#
#   column 1      route 1 — the induced common, ONE partition on ONE plane
#   columns 2..   route 2 — the PCLF, one partition PER NODE, drawn separately because merged they
#                 would superimpose |S| partitions over the same plane and show neither
#   row 1         primal De Bruijn k=1 (COMPLETE): the node records the mode just PLAYED
#   row 2         dual De Bruijn k=1 (CO-COMPLETE): the node commits to the mode played NEXT
#
# The two rows are the campaign's two mechanisms, visible rather than tabulated:
#
#   primal   complete ⟹ the induced common is min_i V_i, whose sublevel set is a non-convex UNION.
#            Column 1 inherits that fragmentation in every cell; columns 2.. keep convex pieces.
#            Route 2 builds MORE cells here and is still faster — cell count is the wrong proxy.
#
#   dual     co-complete ⟹ the induced common is max_i V_i, a convex INTERSECTION — column 1 is one
#            clean piece, so there is no fragmentation for route 2 to avoid. The only gain left is
#            that each node refines under ONE mode instead of two, and whether that pays turns on
#            the per-node reduction `g` against |S| = 2 (plan.md §3a).
#
# ── ON GLB, THE DOMAIN IS THEIRS ──────────────────────────────────────────────────────────────────
#
# Gol, Ding, Lazar & Belta fix their own working set and terminal set in Examples 3.1/4.1:
#
#     X = {x : ‖Lx‖_∞ ≤ Γ_X},  Γ_X = 10      area 308.0
#     D = {x : ‖Lx‖_∞ ≤ Γ_D},  Γ_D = 5.063   area  79.0      (L = `gol_lazar_belta_L()`)
#
# and 11 slices between them at their rate ρ = 0.94. Those are used here, NOT the ±5.9 box that
# `gol_lazar_belta_problem` carries: that box has area 139 and is a strict SUBSET of their X, so a
# quotient built on it tiles well under half their domain and the three observation regions fall
# largely outside it.
#
# THE RUNG COUNT IS DERIVED, NOT CHOSEN, AND IT DIFFERS BETWEEN THE ROWS. The ladder is geometric at
# the certificate's own rate, so the terminal level is pinned to ΓX·γ^(k−1) and `k` is the only knob.
# `k` is taken from the construction's own stopping rule — the first level whose terminal set no
# longer meets R1, R2 or R3 — which is also how their Γ_D = 5.063 is fixed. See
# `rungs_clearing_regions` for why the two rules tried before it were both wrong. It comes out at
# 7 rungs on the primal against 13 on the dual, because the dual contracts more slowly (γ = 0.924
# against 0.867) and needs more levels to reach a region-free terminal set.
#
# WHY ΓX IS COMPUTED PER ROW AND IMPOSED ON BOTH ARMS. Each arm would otherwise derive its outer
# level from its own certificate, and the induced common's gauge is not the PCLF pieces' gauge — the
# two arms would tile different regions and the panels would differ by coverage rather than by
# construction. Route 1's own natural cover of their X is computed once per orientation and imposed
# on route 2: the baseline keeps the geometry it would have chosen, and route 2 is made to conform.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Printf
using Spot
import Statistics

gr()

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)

# Paired rounds per row. Three is the minimum that gives a median; raise it with SAMPLES=n in the
# environment when a row's spreads overlap.
const SAMPLES = parse(Int, get(ENV, "SAMPLES", "3"))

const WHICH = isempty(ARGS) ? "glb" : lowercase(ARGS[1])
WHICH in ("glb", "a") ||
    error("first argument must be \"glb\" or \"A\", got $(repr(ARGS[1]))")

# Their published levels, and the sublevel sets of their own certificate at those levels.
const GLB_GAMMA_X = 10.0
const GLB_GAMMA_D = 5.063
glb_their_sublevel(γ) =
    let L = gol_lazar_belta_L()
        LazySets.HPolytope(vcat(L, -L), fill(γ, 2 * size(L, 1)))
    end

"The problem, its working set, how a PCLF is fitted on it, and how deep to go."
function example(which)
    if which == "glb"
        (; f, R1, R2, R3) = gol_lazar_belta_problem()
        X = glb_their_sublevel(GLB_GAMMA_X)
        return (;
            name = "Gol–Lazar–Belta Example 3.1 (their X and D)",
            slug = "C_gol_lazar_belta",
            problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3]),
            X = X,
            f = f,
            pclf_of = (graph, nodes) -> PCLF.compute_polyhedral_pieces_pclf(
                f,
                graph,
                OPTIMIZER,
                PCLF.conic_partitions_dict_2d(2, nodes);
                MLF = true,
            ),
            rungs = nothing,
            # Their own formula and their own initial point (`a` in the paper). `D` is not listed
            # here: each arm's terminal set is its own certificate's innermost sublevel set, so it is
            # filled in per build.
            spec = (;
                φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))",
                regions = Dict{Symbol, Any}(:R1 => R1, :R2 => R2, :R3 => R3),
                ap_to_obs = Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
                x0 = SVector(-4.0, -7.0),
            ),
        )
    end
    (; f, X, R1, R2) = two_mode_problem()
    return (;
        name = "canonical two-mode example",
        slug = "A_identical_pieces",
        problem = PR.BisimulationQuotientProblem(f, X, [R1, R2]),
        X = X,
        f = f,
        spec = (;
            φ = ltl"(!R2 U R1)",
            regions = Dict{Symbol, Any}(:R1 => R1, :R2 => R2),
            ap_to_obs = Dict(:D => -1, :R1 => 1, :R2 => 2),
            # A point, not the degenerate box `1_comparison.jl` uses: that script has its own runner,
            # while `synthesize_cosafe_ltl` indexes `x0` directly.
            x0 = SVector(1.0, 1.0),
        ),
        pclf_of = (graph, nodes) -> PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
            f,
            graph,
            OPTIMIZER;
            Gmats = rotation_templates(nodes; θ = π / 6, mode = :rotation),
            MLF = true,
            verbose = false,
        ),
        # Derived like GLB's: the canonical example has no published terminal set, but the rule that
        # D must clear the observation regions is the construction's, not GLB's, so it applies here
        # too — the fixed 7 this script used before had never been checked against it.
        rungs = nothing,
    )
end

const EX = example(WHICH)
# An explicit rung count on the command line overrides the derivation, for quick looks.
const FORCED_RUNGS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : EX.rungs

@printf("%s — the four quotients\n\n", EX.name)

common_piece(cert) = cert.pieces[:clf]
sublevel(cert, τ) = PCLF.get_sublevel_set(common_piece(cert), τ; atol = 1e-6)

"""
The first rung count whose terminal set clears every observation region, over all of `certs`.

This is the construction's OWN stopping rule (`all_nodes_clear_regions`), and it is how their Γ_D is
fixed: 5.063 is the first level of their certificate that no longer touches R1, R2 or R3. A terminal
set that still meets a region is not merely coarse — every point of D is labelled with the terminal
observation, so an overlap mislabels those points and the controllable set stops being comparable to
theirs.

Two earlier rules are NOT used and both were wrong:
  * a fixed rung count, which descends further on the faster-contracting graph;
  * `⌈log(Γ_X/Γ_D)/log(1/γ)⌉`, their band DEPTH at our rate — it reproduces the ratio but not the
    property, and it left D overlapping the regions on both orientations (5 rungs where 7 were
    needed on the primal, 9 where 13 were needed on the dual).

The maximum over the arms is taken so that both tile the same ladder: on a co-complete graph the
induced common is the INTERSECTION of the pieces, so it can clear a region at a level where an
individual node's piece does not.
"""
function rungs_clearing_regions(certs, ΓX, γ, regions; kmax = 60, tol = 1e-2)
    ks = map(certs) do cert
        k = findfirst(
            k -> PCQ.all_nodes_clear_regions(cert, ΓX * γ^(k - 1), regions; tol = tol),
            1:kmax,
        )
        isnothing(k) && error("no region-free terminal level within $kmax rungs")
        return k
    end
    return maximum(ks)
end

"""
Check the terminal set that was ACTUALLY produced against every observation region.

`rungs_clearing_regions` selects the rung count from the *level*, which is the construction's own
criterion; this re-checks the set it returned. A `D` meeting a region would label every point of the
overlap with the terminal observation instead of its own, so the controllable sets would stop being
comparable to the reference construction's — it is the one property worth asserting rather than
assuming.
"""
function assert_D_clears_regions(D, regions, label)
    parts = D isa LazySets.UnionSetArray ? D.array : [UT._as_hpolytope(D)]
    for (nm, R) in regions
        RH = UT._as_hpolytope(R)
        for P in parts
            LazySets.isdisjoint(P, RH) || error(
                "[$label] terminal set D meets observation region $nm — the labelling is ambiguous " *
                "and the run is void. Increase the rung count.",
            )
        end
    end
    @printf(
        "    D clears %s (%d part(s))\n",
        join(sort!(string.(collect(keys(regions)))), ", "),
        length(parts)
    )
    return nothing
end

"Solve `spec` on `q`, warmed up then timed over `SAMPLES` rounds."
function timed_solve(f_sys, q, D, spec, label)
    run() = synthesize_cosafe_ltl(
        f_sys,
        q,
        Dionysos.spot_stepper(spec.φ),
        merge(Dict{Symbol, Any}(:D => D), spec.regions),
        spec.ap_to_obs,
        spec.x0;
        print_level = 0,
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
        "    %-28s %8.3f s   %d cells won\n",
        label,
        Statistics.median(ts),
        length(res.controllable_set)
    )
    return (; res, ts)
end

function cell_panel(q, ttl; node = nothing)
    p = plot(; aspect_ratio = :equal, legend = false, title = ttl, titlefontsize = 8)
    opts = (;
        what = :states,
        by = :state,
        show_contours = true,
        linewidth = 0.3,
        fillalpha = 0.9,
        merge_series = false,
    )
    if node === nothing
        plot!(p, q; opts...)
    else
        plot!(p, q; node = node, opts...)
    end
    return p
end

# `PCQ.cell_complexities` returns two unkeyed vectors over the whole quotient, so it cannot answer
# "worst cell of THIS node" — the per-node figure the dual row is about. Read the cells directly.
"Cells, worst-cell parts and worst-cell facets, over the whole quotient or over one node."
function shape(q; node = nothing)
    cells = [v for v in values(q.states) if node === nothing || v.node == node]
    return (;
        n = length(cells),
        parts = maximum(PCQ.num_parts(c.set) for c in cells),
        faces = maximum(PCQ.num_faces(c.set) for c in cells),
    )
end

rows = []
for dual in (false, true)
    label = dual ? "dual" : "primal"
    println("── $label De Bruijn k=1 ", "─"^50)

    graph = PCLF.generate_DeBruijn_edges(2, 1; dual = dual)
    nodes = sort(collect(graph.verts); by = string)
    pclf = EX.pclf_of(graph, nodes)
    common = PCLF.build_common_lyapunov(pclf)
    γ = pclf.JSRapprox

    # The part count is the structural fact the row illustrates: a UNION on the complete graph, a
    # single convex INTERSECTION on the co-complete one.
    common_parts = let S = sublevel(common, 1.0)
        S isa LazySets.UnionSetArray ? length(S.array) : 1
    end

    # `compute_tau_X` is typed on HPolytope; the canonical example declares its working set as a
    # Hyperrectangle, so convert rather than widen a library signature for a script.
    ΓX = PCQ.compute_tau_X(common, UT._as_hpolytope(EX.X))
    nb_levels = if isnothing(FORCED_RUNGS)
        rungs_clearing_regions(
            (common, pclf),
            ΓX,
            γ,
            [UT._as_hpolytope(R) for R in EX.problem.observation_regions],
        )
    else
        FORCED_RUNGS
    end

    @printf(
        "  certified rate %.6f   ΓX = %.4f   %d rungs   induced common: %d part(s)\n",
        γ,
        ΓX,
        nb_levels,
        common_parts
    )

    # Returns `(; quotient, D)`: each arm's terminal set is its OWN certificate's innermost sublevel
    # set, so the specification's `D` atom is per-arm and cannot be hoisted out of here.
    build(cert) = build_quotient(
        EX.problem,
        cert;
        atol = 1e-3,
        nb_levels = nb_levels,
        ΓX = ΓX,
        max_slices = nb_levels,
        print_level = 0,
    )

    # ── TIMING, before a single figure is drawn ───────────────────────────────────────────────────
    # Determinising is INSIDE route 1's clock: it is precisely the step route 2 skips. It is also
    # reported on its own, because "determinising is slow" and "what determinising returns is slow to
    # build on" are different claims and only the second holds.
    route1() = build(PCLF.build_common_lyapunov(pclf))
    route2() = build(pclf)

    println(
        "  warming up (the first call compiles the pipeline; timing it would measure Julia) ...",
    )
    route1()
    route2()

    println("  timing $SAMPLES paired rounds ...")
    t_det, t1, t2 = Float64[], Float64[], Float64[]
    # The last round's quotients are kept for the figures and the shape table. Rebuilding them
    # afterwards would be a third full build per arm — on the deep ladders this problem needs, that
    # is tens of minutes spent reproducing something already in hand.
    b1, b2 = nothing, nothing
    for k in 1:SAMPLES
        # The two arms run back-to-back inside one round, so the ratio is paired: a machine that
        # slows down mid-run moves both arms of that round together instead of one arm of the pool.
        GC.gc()
        a = time_ns()
        cert = PCLF.build_common_lyapunov(pclf)
        b = time_ns()
        b1 = build(cert)
        c = time_ns()
        push!(t_det, (b - a) / 1e9)
        push!(t1, (c - a) / 1e9)

        GC.gc()
        d = time_ns()
        b2 = route2()
        push!(t2, (time_ns() - d) / 1e9)

        @printf(
            "    round %d: route 1 %7.3f s   route 2 %7.3f s   ratio %5.2f\n",
            k,
            t1[k],
            t2[k],
            t1[k] / t2[k]
        )
    end

    q1, q2 = b1.quotient, b2.quotient
    s1, s2 = shape(q1), shape(q2)
    per_node = [(nd, shape(q2; node = nd)) for nd in quotient_nodes(q2)]

    med(v) = Statistics.median(v)
    @printf(
        "  route 1 %7.3f s (min %.3f)  of which determinise %.3f s (%.1f %%)\n",
        med(t1),
        minimum(t1),
        med(t_det),
        100 * med(t_det) / med(t1)
    )
    @printf("  route 2 %7.3f s (min %.3f)\n", med(t2), minimum(t2))
    @printf(
        "  SPEEDUP  paired median %.2fx   minimum-based %.2fx\n",
        med(t1 ./ t2),
        minimum(t1) / minimum(t2)
    )

    @printf(
        "  route 1  %5d cells   worst cell %3d parts / %4d facets\n",
        s1.n,
        s1.parts,
        s1.faces
    )
    @printf(
        "  route 2  %5d cells   worst cell %3d parts / %4d facets   (%s)\n\n",
        s2.n,
        s2.parts,
        s2.faces,
        join([@sprintf("node %s: %d", nd, s.n) for (nd, s) in per_node], ", ")
    )

    # ── THE SPECIFICATION, on both quotients and under both quantifiers ──────────────────────────
    # Synthesis (∃): the modes are the CONTROLLER's — does SOME switching signal satisfy φ?
    # Verification (∀): the modes are the ENVIRONMENT's — does EVERY signal satisfy it?
    # Both are solved because the benchmark poses both, and because a cheaper quotient is only worth
    # having if it answers both questions the same way: the two arms' certified sets are the
    # soundness check on every timing above.
    f_forall = ST.with_switching(EX.f, HybridSystems.AutonomousSwitching())
    solves = Dict{String, Any}()
    for (arm, built) in (("route 1", b1), ("route 2", b2))
        println("  $arm — terminal set and specification")
        assert_D_clears_regions(built.D, EX.spec.regions, "$label / $arm")
        solves["$arm ∃"] =
            timed_solve(EX.f, built.quotient, built.D, EX.spec, "synthesis (∃)")
        solves["$arm ∀"] =
            timed_solve(f_forall, built.quotient, built.D, EX.spec, "verification (∀)")
    end
    @printf(
        "  SOLVE SPEEDUP  ∃ %.2fx   ∀ %.2fx\n\n",
        med(solves["route 1 ∃"].ts ./ solves["route 2 ∃"].ts),
        med(solves["route 1 ∀"].ts ./ solves["route 2 ∀"].ts)
    )

    panels = [
        cell_panel(
            q1,
            "ROUTE 1 · $label, $nb_levels rungs\ninduced common, $common_parts part$(common_parts == 1 ? "" : "s") — $(s1.n) cells, worst $(s1.parts) parts",
        ),
    ]
    for (nd, s) in per_node
        push!(
            panels,
            cell_panel(
                q2,
                "ROUTE 2 · $label — node $nd\n$(s.n) cells, worst $(s.parts) parts";
                node = nd,
            ),
        )
    end
    push!(
        rows,
        (;
            label,
            panels,
            s1,
            s2,
            per_node,
            common_parts,
            ΓX,
            nb_levels,
            t1,
            t2,
            t_det,
            solves,
        ),
    )
end

# One window for all six panels: on separate windows a small covered region and a large one look
# alike, and the two rows really do cover slightly differently.
all_panels = reduce(vcat, [r.panels for r in rows])
lims = align_panels!(all_panels)

ncols = maximum(length(r.panels) for r in rows)
blank() = plot(; framestyle = :none, legend = false)
grid = reduce(
    vcat,
    [vcat(r.panels, [blank() for _ in 1:(ncols - length(r.panels))]) for r in rows],
)

w, h = panel_row_size(lims, ncols)
fig = plot(
    grid...;
    layout = (length(rows), ncols),
    size = (w, length(rows) * h + 34),
    plot_title = "$(EX.name) — route 1 (left) against route 2 (right); primal above, dual below",
    plot_titlefontsize = 10,
)

figdir = joinpath(@__DIR__, EX.slug)
mkpath(figdir)
out = joinpath(figdir, "primal_vs_dual_quotients.png")
savefig(fig, out)
println("wrote $out")

med(v) = Statistics.median(v)

println("\n", "="^92)
@printf(
    "%-8s %6s %7s %8s %8s %9s %6s %9s %9s %8s\n",
    "graph",
    "rungs",
    "common",
    "route1",
    "route2",
    "worst 1/2",
    "g",
    "t route1",
    "t route2",
    "SPEEDUP"
)
println("-"^92)
for r in rows
    per = r.s2.n / length(r.per_node)
    @printf(
        "%-8s %6d %5d p %8d %8d %5d/%-3d %6.2f %8.3fs %8.3fs %7.2fx\n",
        r.label,
        r.nb_levels,
        r.common_parts,
        r.s1.n,
        r.s2.n,
        r.s1.parts,
        r.s2.parts,
        r.s1.n / per,
        med(r.t1),
        med(r.t2),
        med(r.t1 ./ r.t2)
    )
end
println("="^92)

@printf(
    "\n%-8s %-8s %10s %10s %10s %10s %10s %10s\n",
    "graph",
    "arm",
    "build",
    "∃ solve",
    "∃ cells",
    "∀ solve",
    "∀ cells",
    "build+∃+∀"
)
println("-"^92)
for r in rows
    for (arm, t) in (("route 1", r.t1), ("route 2", r.t2))
        e, a = r.solves["$arm ∃"], r.solves["$arm ∀"]
        @printf(
            "%-8s %-8s %9.3fs %9.3fs %10d %9.3fs %10d %9.3fs\n",
            r.label,
            arm,
            med(t),
            med(e.ts),
            length(e.res.controllable_set),
            med(a.ts),
            length(a.res.controllable_set),
            med(t) + med(e.ts) + med(a.ts)
        )
    end
end
println("="^92)
println(
    """
THE ∃/∀ CELL COUNTS ARE NOT COMPARABLE ACROSS ARMS. The two arms partition the same region
differently, so 95 winning cells and 188 winning cells can certify the same area — the count measures
the partition, not the result. The soundness check is on the certified REGION, by volume; it lives in
`1_comparison.jl` and is reported in plan.md §4c, and it is off here because CDDLib over tens of
thousands of semilinear cells costs more than the builds being timed.

What the counts here are good for is a sanity signal within one arm: ∀ must win no more than ∃.

Each arm solves against its OWN `D` — the innermost sublevel set of its own certificate — and every
`D` was asserted disjoint from every observation region before solving.""",
)

println(
    """
SPEEDUP is the PAIRED ratio: the median of route1/route2 taken round by round, not the ratio of the
two medians — a machine that slows down mid-run then moves both arms of a round together. Its
minimum-based twin is printed per row above; read the two together, and treat a row whose per-round
ratios straddle 1 as undecided rather than as a win.

`g` is route 1's cells against route 2's cells PER NODE. On the DUAL row route 2 pays |S| = 2 copies
of the pipeline and `g` is the only saving, so `g > 2` is break-even. On the PRIMAL row `g < 1` is
expected and route 2 wins anyway, on cell SIMPLICITY: read `worst 1/2` (worst-cell parts, route 1
against route 2), not the cell counts.""",
)
