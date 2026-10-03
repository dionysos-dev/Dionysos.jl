# The full comparison in one table: cells, transitions, per-cell complexity, volume, and time.
#
# Cell COUNT is only half the cost. A quotient is expensive in three ways and they do not have to move
# together, so all three are reported side by side:
#
#   how many cells        what the count columns measure
#   how complex each is   parts per cell and facets per cell -- the refinement primitives are
#                         superlinear in facet count, so a cell with 100 facets is far worse than
#                         twenty-five with 4. The DISTRIBUTION matters, not the mean: a long tail is
#                         what actually hurts, which is why the histograms exist.
#   how many transitions  what synthesis walks. A co-safe solve is a fixed point over transitions, not
#                         over cells, so this is the column that predicts downstream cost.
#
# Both cases, because they answer different questions:
#   WITH regions     what a user runs, and the only case where an LTL formula can be posed at all
#                    (the formula names the regions). Timed end to end, including the solve.
#   WITHOUT regions  the pure graph-and-matrices comparison. The ladder must be pinned there (`ΓX` and
#                    `nb_levels`) or the arms tile different sets -- see plan.md.
#
# Volume is carried alongside every count because the two routes DO NOT tile the same region: route 1
# covers the intersection of the pieces' sublevel sets and route 2 their union, so raw counts compare
# a fine partition of a small set against a coarse partition of a big one. Density is the comparable
# quantity and it is reported as such.

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
"""
Two instances of the SAME system, differing only in how small the terminal set `D` is.

The ladder descends by a factor `γ` per rung until it clears the regions, so making `D` smaller means
adding rungs, and each rung multiplies the partition. This is the scaling knob, and the reason to turn
it is that a result measured at one problem size is a data point, not a trend: if route 2's advantage
is a property of the mechanism it must survive a quotient an order of magnitude larger, and if it is a
property of this particular small instance it will not.

- `large_D` is the ladder the regions themselves produce -- 5 rungs, `tauD/tauX = 0.265`. The
  with-regions arm takes `nothing`, which lets the construction stop where it naturally stops; the
  region-free arm is pinned to the 5 rungs that produces, so "with" and "without" differ in ONE
  variable rather than two.
- `small_D` forces 8 rungs in both arms, `tauD/tauX = 0.098`, roughly fifteen times the cells. Fewer
  timing samples because route 1's primal-graph build runs into the minutes there; two still give a
  spread, and the effects being separated are 2x-6x.
"""
const ALL_INSTANCES = (
    (; name = "large_D", with_levels = nothing, bare_levels = 5, samples = 9),
    (; name = "mid_D", with_levels = 7, bare_levels = 7, samples = 5),
    (; name = "small_D", with_levels = 8, bare_levels = 8, samples = 3),
)

"""
Which instances this run measures.

`mid_D` at 7 rungs is the scaling point: each rung multiplies the partition by roughly four, so it
carries 15–20× the cells of `large_D` while staying minutes rather than hours. `small_D` at 8 rungs is
another 4× on top and ran to well over an hour on route 1's primal arm, so it is kept defined and left
out — one rung is the difference between a run you wait for and a run you schedule.

Fewer timing samples above `large_D`, because two still give a spread and the effects being separated
are 2×–8×.
"""
const INSTANCES = filter(i -> i.name in ("large_D", "mid_D"), ALL_INSTANCES)

"""
Whether to compute polytopal volumes (`CDDLib`), which are the slowest thing in this script.

**Off.** Exactness is confirmed *visually* instead, from `fig_certified_*`: both arms build a
bisimulation of the same system for the same formula, so they certify the same region, and the figure
shows the same green set over two very different partitions. At `large_D` the numbers agreed to 2.27 %
(dual) and 1.22 % (primal, normalised by coverage), so there is nothing left for the arithmetic to add.

Two things are given up with it, stated so they are not silently lost:

- the **density** column (cells per unit volume), which is the only count directly comparable between
  arms that tile different regions;
- the **detector** for the arms tiling different regions at all. It matters: at `small_D` with regions
  the covered volumes came out 29.54 against 38.45, 23 % apart, where at `large_D` they agreed to 2-4 %.
  The with-regions case is the one where `ΓX` is NOT pinned, so at 8 rungs the two certificates' gauges
  diverge and their cell counts stop being like-for-like. **At `small_D` the comparable measurement is
  therefore the region-free one, where `ΓX` is pinned.**
"""
const MEASURE_VOLUME = false

"""
Whether this run draws figures, set from the command line: `julia comparison.jl --figures`.

**Timing and plotting must not share a run.** GR allocates heavily and leaves the heap in a different
state, and the figures for one configuration are drawn between the timed builds of the next, so a
plotting run measures the plotting as much as the construction. That is a plausible contributor to the
same configuration reporting 1.55x and 3.11x on consecutive runs.

So there are two passes over exactly the same configurations:

    julia comparison.jl              timing only  — quote the speed-ups from this
    julia comparison.jl --figures    figures only — one sample, no timing claims

Everything structural (cells, transitions, facets, parts) is deterministic and identical in both, so
the two passes can be cross-checked against each other line by line.
"""
const DRAW_FIGURES = "--figures" in ARGS

"Timing samples: the full count for a timing pass, one for a figure pass."
samples_for(inst) = DRAW_FIGURES ? 1 : inst.samples

"""
The graph families every example is run on: both orientations, at De Bruijn orders 1 and 2.

Order 2 is not more of the same. §2c of `plan.md` predicts that route 1's cost tracks the number of
**observer states** rather than its cell count, and the subset construction is exponential in `|S|` in
the worst case — so a larger graph is the one axis on which that account can fail rather than merely
accumulate confirmations. Order 2 quadruples the node count, which also quadruples what route 2 pays
in lifted copies on a complete graph: if the advantage were an artefact of `|S| = 2` it should visibly
erode here.
"""
const GRAPHS = (
    ("dual De Bruijn k=1", true, 1),
    ("primal De Bruijn k=1", false, 1),
    ("dual De Bruijn k=2", true, 2),
    ("primal De Bruijn k=2", false, 2),
)

"""
The two examples the paper uses, and why there must be two.

Route 1's induced common is the **intersection** of the pieces on a co-complete graph — convex
however different they are — and their **union** on a complete one, which fragments as they separate.
So piece diversity costs route 2 on the dual graph and costs route 1 on the primal, and **no single
example can exhibit both regimes**. The shear screen recorded in plan.md §7 established this and
selected `B`.

- **`A_identical_pieces`** — `two_mode_problem` under a shared 4-facet template. The two node pieces
  come out as the *same set* (measured gap 0.0000), so both routes carry identical geometry and the
  entire difference is the graph. A functional gain is structurally impossible on a co-complete graph,
  so nothing can be confounded with the algorithmic one: this is the **control**.
- **`B_diverse_pieces`** — the observer-study pair under a shared conic partition rich enough for the
  solver to *discover* node-dependent shapes instead of having them imposed (gap 0.46, induced common
  of 13 parts). It carries the **same working set, the same two regions and the same formula as A**, so
  the two examples differ in exactly one respect: whether the certificate's node pieces coincide.

Both graphs are run for both examples. On B the dual arm is the one the open questions single out —
route 1's common there is the *intersection*, convex however different the pieces are, so diversity is
expected to cost route 2. Running it is how that expectation gets tested rather than assumed.

One thing to watch on B: its rate is 0.902 against A's 0.718, so its ladder descends far more slowly
and the region-bearing case needs many more rungs to clear the regions. `max_slices` caps it, which
means the with-regions ladder may stop *before* clearing them — so on B the region-free rows remain
the like-for-like ones, and the with-regions rows are the end-to-end user's-eye measurement.
"""
const EXAMPLES = (
    (;
        name = "A_identical_pieces",
        graphs = GRAPHS,
        setup = function ()
            p = two_mode_problem()
            return (;
                p.f,
                p.X,
                regions = [p.R1, p.R2],
                region_names = [:R1, :R2],
                formula = ltl"(!R2 U R1)",
                x0 = LazySets.Hyperrectangle(; low = [1.0, 1.0], high = [1.0, 1.0]),
                D = UT.semilinear_set([UT._as_hpolytope(p.X)]),
            )
        end,
        pclf_of = function (f, graph)
            nodes = sort(collect(graph.verts); by = string)
            return PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
                f,
                graph,
                OPTIMIZER;
                Gmats = rotation_templates(nodes; θ = π / 6, mode = :rotation),
                MLF = true,
                verbose = false,
            )
        end,
    ),
    (;
        name = "B_diverse_pieces",
        graphs = GRAPHS,
        setup = function ()
            A1 = (1.0 / 10.0) * [1.5519 0.4474; 7.6412 7.4716]
            A2 = (1.0 / 10.0) * [0.4750 9.1755; 1.8955 0.1850]
            f = ST.with_switching(
                HybridSystems.discreteswitchedsystem([A1, A2]),
                HybridSystems.ControlledSwitching(),
            )
            # The SAME working set and the SAME two regions as example A, deliberately. The two
            # examples are meant to differ in exactly one respect -- whether the certificate's node
            # pieces coincide -- and giving B its own geometry would add a second difference and make
            # every cross-example reading ambiguous.
            X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
            return (;
                f,
                X,
                regions = [
                    LazySets.Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5]),
                    LazySets.Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5]),
                ],
                region_names = [:R1, :R2],
                formula = ltl"(!R2 U R1)",
                x0 = LazySets.Hyperrectangle(; low = [1.0, 1.0], high = [1.0, 1.0]),
                D = UT.semilinear_set([UT._as_hpolytope(X)]),
            )
        end,
        pclf_of = function (f, graph)
            nodes = sort(collect(graph.verts); by = string)
            return PCLF.compute_polyhedral_pieces_pclf(
                f,
                graph,
                OPTIMIZER,
                PCLF.conic_partitions_dict_2d(2, nodes);
                MLF = true,
            )
        end,
    ),
    (;
        # The published benchmark: Gol, Ding, Lazar & Belta (arXiv:1208.5471, Example 3.1) — their
        # dynamics, their working set, their three observation regions, their co-safe formula and
        # their initial point. Only the certificate is ours.
        #
        # Their own certificate is a COMMON Lyapunov function (‖Lx‖_∞, ρ = 0.94) and comparing
        # against it would be a different experiment: a certificate determines its own slice family
        # and terminal set, so two certificates tile different regions and their cell counts are not
        # comparable. The comparison here is the campaign's — one PCLF from an oracle, then route 1
        # against route 2 on it — run on a problem a reader already knows.
        name = "C_gol_lazar_belta",
        graphs = GRAPHS,
        setup = function ()
            p = gol_lazar_belta_problem()
            return (;
                p.f,
                p.X,
                regions = [p.R1, p.R2, p.R3],
                region_names = [:R1, :R2, :R3],
                formula = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))",
                x0 = SVector(-4.0, -7.0),   # point `a` of the paper
                D = UT.semilinear_set([UT._as_hpolytope(p.X)]),
            )
        end,
        pclf_of = function (f, graph)
            nodes = sort(collect(graph.verts); by = string)
            return PCLF.compute_polyhedral_pieces_pclf(
                f,
                graph,
                OPTIMIZER,
                PCLF.conic_partitions_dict_2d(2, nodes);
                MLF = true,
            )
        end,
    ),
)

"One folder per example per instance, so no two configurations can be read off one another."
function figdir(ex, inst)
    d = joinpath(@__DIR__, ex.name, inst.name)
    mkpath(d)
    return d
end

summarise(ts) =
    (; med = Statistics.median(ts), min = minimum(ts), lo = minimum(ts), hi = maximum(ts))

"""
Time two arms by **alternating** them, one sample of each in turn, and report a PAIRED ratio.

Three things make a speed-up quotable rather than indicative, and the first two were learned by
watching the same configuration report 1.55x and 3.11x on consecutive runs of identical code.

**Paired, not pooled.** The ratio is the median of the per-round ratios `t1[k] / t2[k]`, not the ratio
of the two medians. Because the rounds are interleaved, round `k` of each arm meets nearly the same
machine state, so a common-mode disturbance -- thermal throttling, a background process, a clock
change -- divides out of every paired ratio and does not divide out of a ratio of medians. This is the
whole reason to interleave in the first place, and pooling the samples afterwards threw the benefit
away.

**Minimum, not median, for the absolute figures.** The work is deterministic: the true cost is a
constant and every disturbance can only ADD time. The minimum is therefore the least-contaminated
estimator of the cost itself, which is what `BenchmarkTools` reports and for the same reason. The
median is kept alongside because a large gap between them says the machine was busy.

**The spread of the paired ratios is the honest error bar** -- not the spread of either arm's times,
which mostly measures the machine.
"""
function sample_pair(b1, b2, samples)
    b1()
    b2()
    t1, t2 = Float64[], Float64[]
    for _ in 1:samples, (b, ts) in ((b1, t1), (b2, t2))
        GC.gc()
        t0 = time_ns()
        b()
        push!(ts, (time_ns() - t0) / 1e9)
    end
    paired = t1 ./ t2
    pair = (;
        ratio = Statistics.median(paired),
        ratio_min = minimum(t1) / minimum(t2),
        lo = minimum(paired),
        hi = maximum(paired),
        n = samples,
    )
    return summarise(t1), summarise(t2), pair
end

"A paired speed-up and its spread, printed the way it should be quoted."
fmt_pair(p) = @sprintf(
    "%.2fx [%.2f-%.2f] paired, %.2fx on minima (n=%d)",
    p.ratio,
    p.lo,
    p.hi,
    p.ratio_min,
    p.n
)

fmt(t) = @sprintf("%.2f [%.2f-%.2f]", t.med, t.lo, t.hi)

"""
What one route must do, end to end: from the PCLF it is handed to the quotient it produces.

`pre` is the work that comes BEFORE the build -- route 1's determinisation. It lives inside the closure
so the time reported is the whole bill rather than its second half; an earlier version of this script
computed the common outside and silently handed route 1 a free step.
"""
function builder(problem, pclf; nb_levels = nothing, ΓX = nothing, pre = nothing)
    return function ()
        pre === nothing || pre()
        return build_quotient(
            problem,
            pclf;
            atol = 1e-3,
            nb_levels = nb_levels,
            ΓX = ΓX,
            max_levels = 100,
            max_slices = nb_levels === nothing ? 8 : nb_levels,
            print_level = 0,
        ).quotient
    end
end

"Everything worth knowing about one built quotient, given its measured build time."
function stats_of(q, t)
    parts, faces = PCQ.cell_complexities(q)
    st = PCQ.bisimulation_stats(q)
    covered =
        MEASURE_VOLUME ? PCQ.get_volume(q, keys(q.states); backend = CDDLib.Library()) : NaN
    return (;
        q,
        t,
        cells = st[:num_states],
        transitions = st[:num_transitions],
        nodes = st[:num_nodes],
        slices = st[:num_slices],
        deadends = st[:num_deadend_states],
        faces,
        parts,
        sum_faces = sum(faces),
        mean_faces = Statistics.mean(faces),
        med_faces = Statistics.median(faces),
        max_faces = maximum(faces),
        p90_faces = Statistics.quantile(Float64.(faces), 0.9),
        mean_parts = Statistics.mean(parts),
        max_parts = maximum(parts),
        covered,
        density = st[:num_states] / covered,
    )
end

function table(label, r1, r2, pair = nothing)
    println("\n", "="^100)
    println(label)
    println("="^100)
    rows = [
        ("cells", r1.cells, r2.cells, :lower),
        ("transitions", r1.transitions, r2.transitions, :lower),
        ("nodes", r1.nodes, r2.nodes, :none),
        ("slices", r1.slices, r2.slices, :none),
        ("deadend cells", r1.deadends, r2.deadends, :none),
        ("Σ facets", r1.sum_faces, r2.sum_faces, :lower),
        ("mean facets / cell", r1.mean_faces, r2.mean_faces, :lower),
        ("median facets / cell", r1.med_faces, r2.med_faces, :lower),
        ("90th pct facets / cell", r1.p90_faces, r2.p90_faces, :lower),
        ("MAX facets / cell", r1.max_faces, r2.max_faces, :lower),
        ("mean parts / cell", r1.mean_parts, r2.mean_parts, :lower),
        ("max parts / cell", r1.max_parts, r2.max_parts, :none),
        ("covered volume", r1.covered, r2.covered, :none),
        ("CELLS PER UNIT VOLUME", r1.density, r2.density, :lower),
    ]
    @printf("%-26s %-16s %-16s %s\n", "", "route 1", "route 2", "ratio (r1/r2)")
    println("-"^78)
    # Rows whose quantity was not measured (volume, when MEASURE_VOLUME is off) come through as NaN
    # and are dropped rather than printed as "NaN", which reads like a failure.
    for (name, a, b, dir) in filter(r -> !(r[2] isa AbstractFloat && isnan(r[2])), rows)
        ratio = b == 0 ? "-" : @sprintf("%.2f", a / b)
        mark = dir === :lower && b != 0 ? (a / b > 1.05 ? "  <- route 2 better" : "") : ""
        @printf(
            "%-26s %-16s %-16s %s%s\n",
            name,
            a isa Integer ? string(a) : @sprintf("%.2f", a),
            b isa Integer ? string(b) : @sprintf("%.2f", b),
            ratio,
            mark
        )
    end
    @printf("%-26s %-16s %-16s\n", "build time (s)", fmt(r1.t), fmt(r2.t))
    pair === nothing || @printf("%-26s %s\n", "SPEED-UP (paired)", fmt_pair(pair))
    return nothing
end

"""
Facets-per-cell histograms, which is where the count columns hide the real story.

The mean is a poor summary: the refinement primitives cost superlinearly in facet count, so a long
tail of fat cells dominates even when most cells are simple. Log counts, because most cells sit at the
minimum and the tail is what the eye needs to find.
"""
function histogram_figure(label, r1, r2, stem, dir)
    # The data is plotted as log10(facets) on a LINEAR axis, with ticks relabelled back to facet
    # counts. Two separate reasons, and the second is a trap worth recording.
    #
    # Why log at all: the distribution is heavy-tailed -- 4 facets at the mode, over a thousand in the
    # tail -- so linear bins put ~every cell in the first bar and leave the rest of the plot empty,
    # which is what the first version of this figure did.
    #
    # Why not `xscale = :log10`: it poisons Plots' state for every figure drawn AFTERWARDS. Setting it
    # once made the next equal-aspect figure abort with `total_plotarea_vertical > 0mm`, and in a
    # bisect even a second histogram then failed. Transforming the data instead keeps every axis
    # linear, so no figure in this script depends on the order the figures are drawn in.
    lo = log10(max(1, min(minimum(r1.faces), minimum(r2.faces))))
    hi = log10(max(r1.max_faces, r2.max_faces) * 1.05)
    bins = range(lo, hi; length = 36)
    ticks = filter(
        v -> lo <= log10(v) <= hi,
        [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048],
    )
    fig = plot(;
        size = (900, 430),
        legend = :topright,
        xlabel = "facets per cell",
        ylabel = "number of cells",
        xticks = (log10.(ticks), string.(ticks)),
        left_margin = 6Plots.mm,
        bottom_margin = 5Plots.mm,
        title = label,
        titlefontsize = 10,
    )
    histogram!(
        fig,
        log10.(Float64.(r1.faces));
        bins = bins,
        alpha = 0.55,
        linewidth = 0.3,
        label = "route 1 — $(r1.cells) cells, mean $(round(r1.mean_faces; digits = 1)), max $(r1.max_faces)",
    )
    histogram!(
        fig,
        log10.(Float64.(r2.faces));
        bins = bins,
        alpha = 0.55,
        linewidth = 0.3,
        label = "route 2 — $(r2.cells) cells, mean $(round(r2.mean_faces; digits = 1)), max $(r2.max_faces)",
    )
    savefig(fig, joinpath(dir, "fig_facets_$(stem).png"))
    println("  wrote fig_facets_$(stem).png")
    return nothing
end

"""
A co-safe LTL solve over one quotient, as a closure so the two arms can be timed alternately.

This is the downstream payoff and the reason cell count matters at all: the solve is a fixed point on
the product of the quotient with the specification automaton, so it scales with TRANSITIONS rather
than cells. Only available with regions, since the formula names them.
"""
function ltl_runner(sys, q, D)
    # The formula, the regions and the initial point all come from the example, because a benchmark
    # brings its own specification and rewriting it to fit one hard-coded formula would be measuring
    # a different problem than the one the benchmark poses.
    named = Dict{Symbol, Any}(:D => D)
    for (name, R) in zip(sys.region_names, sys.regions)
        named[name] = R
    end
    ap_to_obs = Dict{Symbol, Int}(:D => -1)
    for (i, name) in enumerate(sys.region_names)
        ap_to_obs[name] = i
    end
    return function ()
        p = PR.CoSafeLTLProblem(
            sys.f,
            sys.x0,
            Dionysos.spot_stepper(sys.formula),
            named,
            Dict{Symbol, Any}(ap => MP.INNER for ap in keys(named)),
        )
        o = MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
        MOI.set(o, MOI.RawOptimizerAttribute("concrete_problem"), p)
        MOI.set(o, MOI.RawOptimizerAttribute("bisimulation_quotient"), q)
        MOI.set(o, MOI.RawOptimizerAttribute("ap_to_obs"), ap_to_obs)
        MOI.set(o, MOI.RawOptimizerAttribute("early_stop"), false)
        MOI.set(o, MOI.RawOptimizerAttribute("print_level"), 0)
        MOI.optimize!(o)
        return MOI.get(o, MOI.RawOptimizerAttribute("controllable_set"))
    end
end

"""
The certified region of one solve, as a VOLUME rather than a cell count.

Every arm builds a bisimulation, so each certifies the EXACT set of states satisfying the formula --
the arms must agree as REGIONS however different their partitions are. Cell counts cannot show that:
174 cells and 379 cells can cover the same area. Only the volume can, and if the volumes disagree
beyond the tolerated sources, exactness is lost and every cost number in this folder is void.
"""
function certified(q, won, t)
    vol = if !MEASURE_VOLUME
        NaN
    elseif isempty(won)
        0.0
    else
        PCQ.get_volume(q, won; backend = CDDLib.Library())
    end
    return (; t, won, vol, ncells = length(won))
end

"""
The certificate itself, before any quotient exists: the PCLF's node pieces, then the common induced
from them. Two panels, one figure, per graph.

This is the picture behind §3b, and it is the only figure here that explains *why* rather than showing
*what*. `ObserverCLFPiece` defines the induced common as

    V*(x) = min over observer states S of ( max over nodes i in S of V_i(x) )

so its sublevel set is a **UNION over observer states of INTERSECTIONS over the nodes in each**. The
number of parts is therefore the number of reachable subsets, and the two graphs sit at the extremes:

  dual    node j enables only mode j but leads to EVERY node, so reading a mode never narrows the
          subset -- the observer is absorbing, one state, and the common is the single convex
          INTERSECTION of the two pieces.
  primal  node i records the last mode, so reading mode m pins the subset to {(m)} -- the observer
          separates into singletons, three states, and the common is the non-convex UNION of the two
          pieces.

Which is the counter-intuitive part worth seeing: the determinisation that *succeeds* in separating
states is the one that produces the more complicated certificate, and every cell route 1 later builds
inherits that shape.
"""
function certificate_figure(gname, pclf, common, dir, stem; γ = 1.0, nlevels = 5)
    nodes = sort(collect(pclf.graph.verts); by = string)
    parts_of(S) = S isa LazySets.UnionSetArray ? S.array : [S]
    nfacets(P) = length(LazySets.constraints_list(P))

    # SEVERAL nested levels, at the certificate's OWN rate. The construction descends the ladder
    # τ, ρτ, ρ²τ, … and tiles the annuli between consecutive rungs, so these are the actual slices the
    # quotient is built on rather than an arbitrary choice of contour. Drawn largest first, so each
    # smaller level lands on top and the nesting reads.
    ρ = pclf.JSRapprox
    levels = [γ * ρ^(k - 1) for k in 1:nlevels]
    shade = Plots.palette(:viridis, nlevels)

    # ONE PANEL PER PIECE, not the pieces overlaid. On this system the solver returns nearly the same
    # polytope at every node, so superimposing them renders as a single grey blob that shows neither --
    # and "the pieces are nearly identical" is itself a result worth being able to read off the figure,
    # since it is what rules out certificate GEOMETRY as the cause of the route-1/route-2 gap and
    # leaves graph structure as the only remaining explanation.
    panels = Any[]
    for nd in nodes
        outer = PCLF.get_sublevel_set(pclf.pieces[nd], first(levels))
        fig = plot(;
            aspect_ratio = :equal,
            legend = false,
            title = "node $nd — $(nfacets(outer)) facets per level",
            titlefontsize = 9,
        )
        for (k, lv) in enumerate(levels)
            plot!(
                fig,
                PCLF.get_sublevel_set(pclf.pieces[nd], lv);
                alpha = 1.0,
                linewidth = 1.2,
                linecolor = :black,
                color = shade[k],
            )
        end
        push!(panels, fig)
    end

    # The common at every level too, each level's parts drawn in that level's colour. The part COUNT is
    # the message and it is reported in the title; colouring by part instead would make the nesting
    # unreadable, and the nesting is what shows that the fragmentation repeats at every rung rather
    # than being an artefact of one contour.
    counts = Int[]
    comfig = plot(; aspect_ratio = :equal, legend = false, titlefontsize = 9)
    for (k, lv) in enumerate(levels)
        cparts = parts_of(PCLF.get_sublevel_set(common.pieces[:clf], lv; atol = 1e-6))
        push!(counts, length(cparts))
        for P in cparts
            plot!(
                comfig,
                P;
                alpha = 1.0,
                linewidth = 1.2,
                linecolor = :black,
                color = shade[k],
            )
        end
    end
    npart = maximum(counts)
    plot!(
        comfig;
        title = "the INDUCED COMMON — $(npart) part$(npart == 1 ? "" : "s") per level",
    )
    push!(panels, comfig)

    # How different the pieces actually are, as the largest relative gap between their support
    # functions. This instance is the one every number in plan.md rests on, and the shear screen of
    # plan.md §7 found a family where route 2's advantage appears ONLY where the gap is zero -- so THIS
    # instance sits on that axis decides whether the headline result is in the degenerate regime or
    # not. Cheap, exact for convex sets, and it belongs next to the figure rather than in a
    # throwaway script.
    gap = if length(nodes) == 2
        outer = [PCLF.get_sublevel_set(pclf.pieces[nd], first(levels)) for nd in nodes]
        maximum(let d = [cos(2π * k / 64), sin(2π * k / 64)]
                h1, h2 = LazySets.ρ(d, outer[1]), LazySets.ρ(d, outer[2])
                abs(h1 - h2) / max(abs(h1), abs(h2), eps())
            end for k in 0:63)
    else
        NaN
    end
    @printf(
        "  piece gap (%s): %.4f  — 0 would mean the two node pieces coincide\n",
        gname,
        gap
    )

    lims = align_panels!(panels)
    fig = plot(
        panels...;
        layout = (1, length(panels)),
        size = panel_row_size(lims, length(panels); width = 400, chrome = 115),
        plot_title = "$gname — $nlevels nested sublevel sets (ρ = $(round(ρ; digits = 4)))   ·   the $(length(nodes)) pieces, then the common they induce",
        plot_titlefontsize = 10,
    )
    savefig(fig, joinpath(dir, "fig_certificate_$(stem).png"))
    println("  wrote fig_certificate_$(stem).png   (parts per level: ", counts, ")")
    return nothing
end

"""
The partitions themselves: route 1's single plane, then route 2's `|S|` lifted copies.

Per-STATE colours, not per-slice. Colouring by slice draws the sublevel annuli -- handsome, and it
hides the only thing the figure exists to show, since every cell in a slice takes one colour and the
379-against-174 difference becomes invisible.

Route 2's nodes are drawn separately. Merged they would superimpose `|S|` partitions over one plane
and show neither; separated they ARE the argument -- one finely cut plane against `|S|` coarse ones.
"""
function partition_figure(gname, r1, r2, stem, dir)
    panel(q, title; node = nothing) = begin
        fig = plot(; aspect_ratio = :equal, legend = false, title = title, titlefontsize = 9)
        plot!(
            fig,
            q;
            what = :states,
            by = :state,
            node = node,
            show_contours = true,
            linewidth = 0.3,
            fillalpha = 0.9,
            merge_series = false,
        )
        return fig
    end

    panels = Any[panel(r1.q, "route 1 — induced common\n$(r1.cells) cells, one plane")]
    for nd in quotient_nodes(r2.q)
        n_here = count(s -> s.node == nd, values(r2.q.states))
        push!(panels, panel(r2.q, "route 2 — $gname\nnode $nd, $n_here cells"; node = nd))
    end
    lims = align_panels!(panels)
    fig = plot(
        panels...;
        layout = (1, length(panels)),
        size = panel_row_size(lims, length(panels)),
    )
    savefig(fig, joinpath(dir, "fig_partitions_$(stem).png"))
    println("  wrote fig_partitions_$(stem).png")
    return nothing
end

"""
The certified region of each arm, on one shared window.

The figure that VALIDATES the comparison rather than illustrating it: the same green region over
wildly different partitions. Legitimate sources of disagreement, which must not be mistaken for a
defect: the covered regions differ by 2-4 %, `atol` erodes a sliver at every cut (and route 1 makes
far more cuts on the primal graph), and the two graphs' certificates differ slightly in rate
(0.717572 against 0.715900), so their slice families do not coincide exactly.
"""
function certified_figure(entries, stem, dir)
    panels = map(entries) do (label, q, won)
        fig = plot(; aspect_ratio = :equal, legend = false, title = label, titlefontsize = 8)
        losing = collect(setdiff(Set(keys(q.states)), Set(won)))
        _plot_winning_losing!(
            fig,
            q,
            collect(won),
            losing,
            nothing;
            show_contours = true,
            linewidth = 0.4,
        )
        return fig
    end
    lims = align_panels!(panels)
    fig = plot(
        panels...;
        layout = (1, length(panels)),
        size = panel_row_size(lims, length(panels); width = 370),
    )
    savefig(fig, joinpath(dir, "fig_certified_$(stem).png"))
    println("  wrote fig_certified_$(stem).png")
    return nothing
end

"""
Cells per slice and the growth factor between consecutive slices, per arm.

The totals say route 2 wins; this says *how the win accumulates*, and it is the sequence a closed form
would have to reproduce rather than just its sum. Only meaningful without regions and with the ladder
pinned: the arms then share their levels exactly, so slice `i` means the same annulus in both.

Route 2's counts are divided by `|S|`, because the mechanism is a statement about ONE node's
partition. The last column must exceed `|S|` for route 2 to win overall.
"""
function slice_profile(r1, r2)
    S = r2.nodes
    by1, by2 = Dict(PCQ.states_by_slice(r1.q)), Dict(PCQ.states_by_slice(r2.q))
    n = max(PCQ.num_slices(r1.q), PCQ.num_slices(r2.q))
    c1 = [get(by1, i, 0) for i in 1:n]
    c2 = [get(by2, i, 0) / S for i in 1:n]
    println("\n  per-slice growth (route 2 divided by |S| = $S):")
    @printf(
        "  %-7s %-12s %-10s %-14s %-10s %s\n",
        "slice",
        "r1 cells",
        "growth",
        "r2 per node",
        "growth",
        "r1 / r2-per-node"
    )
    println("  ", "-"^74)
    for i in 1:n
        g1 = i == 1 || c1[i - 1] == 0 ? NaN : c1[i] / c1[i - 1]
        g2 = i == 1 || c2[i - 1] == 0 ? NaN : c2[i] / c2[i - 1]
        ratio = c2[i] == 0 ? NaN : c1[i] / c2[i]
        @printf(
            "  %-7d %-12d %-10s %-14.1f %-10s %s\n",
            i,
            c1[i],
            isnan(g1) ? "-" : @sprintf("%.2f", g1),
            c2[i],
            isnan(g2) ? "-" : @sprintf("%.2f", g2),
            isnan(ratio) ? "-" : @sprintf("%.2fx", ratio)
        )
    end
    return nothing
end

# ---------------------------------------------------------

"""
One oracle graph of one example, measured in both cases.

`dual` and `primal` are run through exactly the same code so the contrast is structural and not a
difference in how each was set up: same system, same template, same tolerances, same ladder depth.
On example A the primal graph is the CONTROL -- its node records the past and still serves every mode,
so its per-node partition should match route 1's and route 2 should lose by roughly |S|.

An example with no observation regions runs the region-free case only: the with-regions case exists to
measure what a user pays on a specification, and without regions there is no specification to pose.
"""
function both_cases(ex, sys, gname, dual, order, inst)
    with_problem = PR.BisimulationQuotientProblem(sys.f, sys.X, sys.regions)
    bare = PR.BisimulationQuotientProblem(sys.f, sys.X, typeof(sys.X)[])
    has_regions = !isempty(sys.regions)

    p2 = ex.pclf_of(sys.f, PCLF.generate_DeBruijn_edges(2, order; dual = dual))
    (isfinite(p2.JSRapprox) && p2.JSRapprox < 1) ||
        (println("\n", gname, ": no certificate"); return nothing)
    p1 = PCLF.build_common_lyapunov(p2)
    isapprox(p1.JSRapprox, p2.JSRapprox; rtol = 1e-4) ||
        (println("\n", gname, ": rates differ, not a cost comparison"); return nothing)

    println("\n\n", "#"^100)
    @printf(
        "# ORACLE PCLF: %s   (rate %.6f, |S| = %d)   INSTANCE: %s\n",
        gname,
        p2.JSRapprox,
        length(p2.graph.verts),
        inst.name
    )
    println("#"^100)

    stem = (startswith(gname, "dual") ? "dual" : "primal") * "_k$(order)"
    dir = figdir(ex, inst)

    # Drawn once per graph, from the certificate alone -- it does not depend on the instance, so it
    # only needs the first one.
    DRAW_FIGURES &&
        inst === first(INSTANCES) &&
        certificate_figure(gname, p2, p1, dir, stem)

    # Route 1's build is timed WITH its determinisation, because that is what a user pays: handed a
    # PCLF, route 1 must produce the common before it can build anything. Excluding it -- which an
    # earlier version of this script did -- silently hands route 1 a free step. On the dual graph it
    # costs nothing (|S| = 2, the observer construction is trivial), but the primal graph's common is a
    # five-part union and is not free.
    p1_of() = PCLF.build_common_lyapunov(p2)

    # The with-regions case exists to measure what a user pays on a specification, so an example
    # that defines no regions has nothing to pose and skips it.
    w1 = w2 = s1 = s2 = wpair = spair = nothing
    if has_regions
        bw1 = builder(with_problem, p1; nb_levels = inst.with_levels, pre = p1_of)
        bw2 = builder(with_problem, p2; nb_levels = inst.with_levels)
        tw1, tw2, wpair = sample_pair(bw1, bw2, samples_for(inst))
        w1, w2 = stats_of(bw1(), tw1), stats_of(bw2(), tw2)
        table("$gname — WITH observation regions [$(inst.name)]", w1, w2, wpair)
        DRAW_FIGURES && histogram_figure(
            "$gname, with regions ($(inst.name))",
            w1,
            w2,
            "$(stem)_with_regions",
            dir,
        )
        DRAW_FIGURES && partition_figure(gname, w1, w2, "$(stem)_with_regions", dir)

        println("\n  co-safe LTL solve of ", sys.formula, ":")
        run1 = ltl_runner(sys, w1.q, sys.D)
        run2 = ltl_runner(sys, w2.q, sys.D)
        # The initial set is one point, and the quotient does not always contain it: at De Bruijn
        # order 2 the certified region shrinks, and example B's dual k=2 quotient misses (1, 1)
        # entirely ("Initial set does not intersect any quotient state"). That is a property of the
        # instance, not a defect, and it must not abort the other configurations — one unsatisfiable
        # specification is a row to report, not a run to lose.
        solvable = try
            run1()
            run2()
            true
        catch err
            @warn "co-safe solve skipped — the quotient does not contain the initial set" gname inst =
                inst.name err = sprint(showerror, err)
            false
        end
        if solvable
            ts1, ts2, spair = sample_pair(run1, run2, samples_for(inst))
            s1, s2 = certified(w1.q, run1(), ts1), certified(w2.q, run2(), ts2)
            for (who, s, w) in (("route 1", s1, w1), ("route 2", s2, w2))
                @printf(
                    "  %-26s %-16s over %-6d transitions   %d cells certified%s\n",
                    who,
                    fmt(s.t),
                    w.transitions,
                    s.ncells,
                    MEASURE_VOLUME ? @sprintf(", volume %.4f", s.vol) : ""
                )
            end
            # THE exactness check, when volumes are measured: both arms build a bisimulation of the same system
            # for the same formula, so they certify the same region -- different partitions, identical answer.
            # TWO numbers, because the raw one is not like-for-like: the arms do not tile the same region (route
            # 1 covers the pieces' sublevel INTERSECTION and route 2 their UNION) and `atol` insets every cut, so
            # the arm making more cuts loses more area. What must agree is the certified FRACTION of what each
            # arm actually tiles; the raw difference is printed beside it so the erosion stays visible.
            # With volumes off this is a VISUAL check instead -- see `MEASURE_VOLUME`.
            if MEASURE_VOLUME
                gap = abs(s1.vol - s2.vol) / max(s1.vol, eps())
                fr1, fr2 = s1.vol / w1.covered, s2.vol / w2.covered
                fgap = abs(fr1 - fr2) / max(fr1, eps())
                @printf(
                    "  %-26s raw %.2f %% apart   |   %.4f vs %.4f of each arm's covered area, %.2f %% apart  %s\n",
                    "certified regions agree?",
                    100 * gap,
                    fr1,
                    fr2,
                    100 * fgap,
                    fgap < 0.05 ? "OK" : "*** DISAGREE ***"
                )
                fgap < 0.05 || @warn(
                    "the two routes certify DIFFERENT regions -- exactness is lost and the cost numbers are void",
                    gname,
                    route1 = s1.vol,
                    route2 = s2.vol,
                    fraction1 = fr1,
                    fraction2 = fr2
                )
            else
                @printf(
                    "  %-26s check the figure: the same green set over two different partitions\n",
                    "certified regions agree?"
                )
            end
            DRAW_FIGURES && certified_figure(
                [
                    (
                        "route 1 — induced common\n$(s1.ncells)/$(w1.cells) cells",
                        w1.q,
                        s1.won,
                    ),
                    ("route 2 — $gname\n$(s2.ncells)/$(w2.cells) cells", w2.q, s2.won),
                ],
                stem,
                dir,
            )
            @printf(
                "  %-26s solve %.2fx, transitions %.2fx\n",
                "ratio (r1/r2)",
                spair.ratio,
                w1.transitions / w2.transitions
            )
        end
    end

    # The ladder is pinned to route 1's own natural tauX so both arms tile the same set -- and to the
    # SAME depth the regions produced, so "without regions" differs from "with" in one variable only.
    nb = inst.bare_levels
    probe = build_quotient(
        bare,
        p1;
        atol = 1e-3,
        nb_levels = nb,
        max_slices = nb,
        print_level = 0,
    )
    ΓX = maximum(MOI.get(probe.optimizer, MOI.RawOptimizerAttribute("Γ")))
    bb1 = builder(bare, p1; nb_levels = nb, ΓX = ΓX, pre = p1_of)
    bb2 = builder(bare, p2; nb_levels = nb, ΓX = ΓX)
    tb1, tb2, bpair = sample_pair(bb1, bb2, samples_for(inst))
    b1, b2 = stats_of(bb1(), tb1), stats_of(bb2(), tb2)
    table(
        "$gname — WITHOUT regions [$(inst.name)] ($nb rungs, ΓX = $(round(ΓX; digits = 4)))",
        b1,
        b2,
    )
    slice_profile(b1, b2)
    DRAW_FIGURES && histogram_figure(
        "$gname, no regions ($(inst.name))",
        b1,
        b2,
        "$(stem)_no_regions",
        dir,
    )
    DRAW_FIGURES && partition_figure(gname, b1, b2, "$(stem)_no_regions", dir)
    return (; ex, gname, stem, inst, w1, w2, b1, b2, s1, s2, wpair, bpair, spair)
end

results = NamedTuple[]
for ex in EXAMPLES
    sys = ex.setup()
    for inst in INSTANCES
        println("\n\n", "#"^100)
        println("# EXAMPLE $(uppercase(ex.name))   ·   INSTANCE $(uppercase(inst.name))")
        println("#"^100)
        for (gname, dual, order) in ex.graphs
            r = both_cases(ex, sys, gname, dual, order, inst)
            r === nothing || push!(results, r)
        end
    end
end

println("\n\n", "="^100)
println("SUMMARY — route 1 / route 2, so above 1 means route 2 is better")
println("="^100)
@printf(
    "%-20s %-10s %-22s %-14s %-9s %-12s %-12s %-10s %s\n",
    "example",
    "instance",
    "oracle graph",
    "case",
    "cells",
    "transitions",
    "max facets",
    "Σ facets",
    "build time"
)
println("-"^130)
for r in results
    # `w1` is `nothing` for an example that defines no regions, so the with-regions row is skipped
    # rather than printed as zeros.
    rows =
        r.w1 === nothing ? (("no regions", r.b1, r.b2, r.bpair),) :
        (("with regions", r.w1, r.w2, r.wpair), ("no regions", r.b1, r.b2, r.bpair))
    for (k, (case, a, b, p)) in enumerate(rows)
        # The speed-up is the PAIRED statistic with its spread, never a ratio of medians — see
        # `sample_pair`. On a figure pass there is one sample, so the spread is degenerate and the
        # column is there only to line up with the timing pass.
        @printf(
            "%-20s %-10s %-22s %-14s %-9.2f %-12.2f %-12.2f %-10.2f %s\n",
            k == 1 ? r.ex.name : "",
            k == 1 ? r.inst.name : "",
            k == 1 ? r.gname : "",
            case,
            a.cells / b.cells,
            a.transitions / b.transitions,
            a.max_faces / b.max_faces,
            a.sum_faces / b.sum_faces,
            @sprintf("%.2f [%.2f-%.2f]", p.ratio, p.lo, p.hi)
        )
    end
end

# The scaling question this second instance exists to answer, stated as a number rather than left to
# the reader: does the advantage hold when the quotient grows by an order of magnitude?
println("\n", "="^100)
println("SCALING — the same ratio at the smallest and largest problem size measured")
println("="^100)
# The two ends of whatever `INSTANCES` holds, rather than two hard-coded names: the point is the
# span between the cheapest and the most expensive configuration actually run.
const SMALLEST, LARGEST = first(INSTANCES).name, last(INSTANCES).name
@printf(
    "%-22s %-14s %-9s %-16s %-14s %-14s %s\n",
    "oracle graph",
    "case",
    "quantity",
    "route 1 cells",
    SMALLEST,
    LARGEST,
    "holds?"
)
println("-"^104)
find(ename, gname, iname) =
    findfirst(r -> r.ex.name == ename && r.gname == gname && r.inst.name == iname, results)
for ex in EXAMPLES,
    (gname, _, _) in ex.graphs,
    (case, sel) in (
        ("with regions", r -> (r.w1, r.w2, r.wpair)),
        ("no regions", r -> (r.b1, r.b2, r.bpair)),
    )

    SMALLEST == LARGEST && continue
    i, j = find(ex.name, gname, SMALLEST), find(ex.name, gname, LARGEST)
    (i === nothing || j === nothing) && continue
    (a1, c1, p1) = sel(results[i])
    (a2, c2, p2) = sel(results[j])
    (a1 === nothing || a2 === nothing) && continue
    for (what, big, small) in
        (("time", p1.ratio, p2.ratio), ("cells", a1.cells / c1.cells, a2.cells / c2.cells))
        @printf(
            "%-22s %-14s %-9s %-16s %-14.2f %-14.2f %s\n",
            what == "time" ? gname : "",
            what == "time" ? case : "",
            what,
            "$(a1.cells) -> $(a2.cells)",
            big,
            small,
            small >= big ? "yes, grows" : (small > 1 ? "yes, shrinks" : "no")
        )
    end
end

"""
The two graphs' arms side by side — the comparison the paired rows structurally cannot show.

Route 1 runs on a DIFFERENT input in each row (the common induced from the dual PCLF is not the one
induced from the primal PCLF), so "how different are those two baselines?" is a real question and not a
re-slicing of the same numbers. The answer is the whole of §3b: they build almost the SAME NUMBER OF
CELLS out of wildly different material. That is the cleanest available demonstration that cell count is
the wrong cost proxy -- here it is held nearly fixed while everything that matters moves by 4x-8x.

Ratios are primal/dual, so above 1 means the primal graph is worse for that arm.
"""
function cross_graph_table(ename, iname, order, case_name, sel)
    i = findfirst(
        r -> r.ex.name == ename && r.stem == "dual_k$(order)" && r.inst.name == iname,
        results,
    )
    j = findfirst(
        r -> r.ex.name == ename && r.stem == "primal_k$(order)" && r.inst.name == iname,
        results,
    )
    (i === nothing || j === nothing) && return nothing
    for (arm, pick) in (
        ("route 1 — the INDUCED COMMON", x -> x[1]),
        ("route 2 — the PCLF directly", x -> x[2]),
    )
        a, b = pick(sel(results[i])), pick(sel(results[j]))
        (a === nothing || b === nothing) && continue
        println("\n", "="^100)
        println("$arm, dual vs primal   [$iname, $case_name]")
        println("="^100)
        @printf(
            "%-26s %-18s %-18s %s\n",
            "",
            "from dual",
            "from primal",
            "ratio (primal/dual)"
        )
        println("-"^82)
        for (name, x, y) in (
            ("cells", a.cells, b.cells),
            ("transitions", a.transitions, b.transitions),
            ("Σ facets", a.sum_faces, b.sum_faces),
            ("mean facets / cell", a.mean_faces, b.mean_faces),
            ("median facets / cell", a.med_faces, b.med_faces),
            ("90th pct facets / cell", a.p90_faces, b.p90_faces),
            ("MAX facets / cell", a.max_faces, b.max_faces),
            ("mean parts / cell", a.mean_parts, b.mean_parts),
            ("max parts / cell", a.max_parts, b.max_parts),
        )
            @printf(
                "%-26s %-18s %-18s %.2f\n",
                name,
                x isa Integer ? string(x) : @sprintf("%.2f", x),
                y isa Integer ? string(y) : @sprintf("%.2f", y),
                y / x
            )
        end
        @printf(
            "%-26s %-18s %-18s %.2f\n",
            "build time (s)",
            fmt(a.t),
            fmt(b.t),
            b.t.min / a.t.min
        )
    end
    return nothing
end

for ex in EXAMPLES,
    inst in INSTANCES,
    (case_name, sel) in (
        ("with regions", r -> (r.w1, r.w2, r.wpair)),
        ("no regions", r -> (r.b1, r.b2, r.bpair)),
    )

    for order in unique(g[3] for g in ex.graphs)
        cross_graph_table(ex.name, inst.name, order, case_name, sel)
    end
end

println(
    """

Reading it. Every column is route 1 divided by route 2, so **above 1 means route 2 is better**.
`build time` includes route 1's determinisation. The primal graph is the CONTROL and should come out at
or below 1 on the count columns -- its node records the past and still serves every mode, so lifting
buys nothing.

One caveat on the count columns, since volumes are not measured here (see `MEASURE_VOLUME`): the two
routes tile slightly different regions, route 1 the pieces' sublevel INTERSECTION and route 2 their
UNION. At `large_D` that difference was measured at 2-4 %, small enough to ignore; with regions at
`small_D` it reached 23 %, so there the region-free rows -- where `ΓX` is pinned and both arms start
from the same tau_X -- are the ones to trust.""",
)
