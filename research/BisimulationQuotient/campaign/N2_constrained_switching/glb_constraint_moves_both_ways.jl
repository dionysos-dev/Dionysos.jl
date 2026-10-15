# Experiment 1: the constraint moves synthesis and verification in OPPOSITE directions.
#
# Run on the **unmodified** Gol-Lazar-Belta example -- their matrices, working set, three regions and
# co-safe formula -- with and without the language "no two consecutive uses of mode 1". Their system
# is stable under arbitrary switching, so their construction works here and E8's numbers are the
# baseline. That is the point: this experiment validates the constrained machinery against a known
# answer rather than asserting a new one.
#
# The prediction, and it is a distinctive one:
#
#   synthesis   (exists)  the CONTROLLER owns the switching. Forbidding 11 removes options from the
#                         controller, so it can do strictly less  ->  winning set SHRINKS.
#   verification (forall) the ADVERSARY owns the switching. Forbidding 11 removes options from the
#                         adversary, so fewer behaviours can violate  ->  verified set GROWS.
#
# Same constraint, same quotient, opposite effects. A constraint that moved both the same way would
# mean the exists/forall distinction is not being handled correctly, so this is a real test.
#
# Baseline from E8 (unconstrained, 2-node De Bruijn, conic order 2): 8 794 of 10 611 cells certified
# for synthesis, 2 962 verified.
#
# Note the constraint here is NOT a dwell time. "No 11" is a maximum run-length on mode 1; a minimum
# dwell time is the opposite and would be the wrong instrument for a system whose single mode is the
# problem. See `dwell_time_example.jl` for the other mechanism.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Spot
using LinearAlgebra

const ORDER = 2

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()

# Clarabel's DEFAULT tolerances, as `examples/gol_lazar_belta_pclf.jl` uses -- deliberately, and it
# matters far more than it looks. Loosening them to 1e-6 (the setting N4 needs at conic order 3,
# where the defaults fail to converge) produced the SAME rate to six digits but a quotient of
# 31 000+ cells against E8's 10 611, because a sloppier solve returns less precise P matrices, so the
# sublevel sets are slightly off and refinement has to work much harder. At order 2 the defaults
# converge fine -- N4 measures 0.867388 against 0.867381 -- so the loosened setting bought nothing
# and cost 3x the cells. Only loosen where order 3 forces it.
const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 1000,
    "verbose" => false,
)

φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"

# Their initial point `a`, exactly as `examples/gol_lazar_belta_pclf.jl` uses it.
const X0 = SVector(-4.0, -7.0)

"""
The certified set under `system`, without building a controller.

`synthesize_cosafe_ltl` also constructs a controller and simulates it, which needs an initial state
that is both covered by the quotient and winning -- neither guaranteed once `max_slices` truncates
the covered region, and a point inside the terminal set is covered but cannot discharge `F(R1)`.
Only the certified set is compared here, so the synthesis stops before the controller.
"""
function certified_set(system, quotient, regions, ap_to_obs; semantics::Symbol = :language)
    # The initial set must be a POINT, as `gol_lazar_belta_pclf.jl` passes. It is not merely an
    # entry point for a trajectory: the co-safe product's initial MEMORY state is fixed by the
    # observation there, and `project_initial_memory_controllable_set` reports the certified set at
    # that memory. Passing the whole working set, which spans all three regions at once, starts the
    # monitor from the wrong state -- and silently: an earlier version of this script did exactly
    # that and reported 442 of 10 611 cells certified where E8 reports 8 794, on an identical
    # quotient. The 20x gap was the only sign anything was wrong.
    initial_set = LazySets.Hyperrectangle(; low = [X0[1], X0[2]], high = [X0[1], X0[2]])
    problem = PR.CoSafeLTLProblem(
        system,
        initial_set,
        Dionysos.spot_stepper(φ),
        regions,
        Dict{Symbol, Any}(ap => MP.INNER for ap in keys(regions)),
    )
    optimizer = MOI.instantiate(PCQ.OptimizerCoSafeLTLOnQuotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), problem)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"), quotient)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("ap_to_obs"), ap_to_obs)
    # `early_stop = false` is what makes this the certified set of the WHOLE domain rather than of
    # the initial set alone. With `true` -- the optimizer's default, and the setting this script
    # omitted at first -- the product is seeded only from states meeting the initial set, so the
    # fixed point explores outward from x0 and stops: 8 cells where E8 reports 8 794, on an
    # identical quotient. `synthesize_cosafe_ltl` passes `false` for exactly this reason.
    MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
    # Scenario B: the plant itself cannot play `11` -- a scheduler, an interlock. The graph's
    # language IS the set of behaviours, so a node that enables no mode-1 edge is describing a move
    # the adversary cannot make, not one the abstraction mislaid. Under the default `:arbitrary`
    # every such gap is charged to the adversary as a loss, which is correct when the graph is only
    # a proof device and catastrophic here: it sank all 446 node-1 cells and, since `∀` fails a
    # state when any successor fails, emptied the answer backwards from there.
    MOI.set(optimizer, MOI.RawOptimizerAttribute("switching_semantics"), semantics)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 0)
    MOI.optimize!(optimizer)
    return (
        set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set")),
        completions = MOI.get(optimizer, MOI.RawOptimizerAttribute("num_completions")),
    )
end

"""
Build the quotient under `graph`, then answer their formula both existentially and universally.

The universal answer comes from re-solving with the switching declared autonomous, which is how the
folder asks the `forall` question -- the modes become the environment's rather than the controller's.
"""
function both_semantics(name, graph, cache_file)
    # The certificate itself costs seconds, so it is recomputed even on a cache hit: the rate is
    # reported in every table here and a cached number that no longer matches the code would be
    # worse than the seconds saved.
    pclf = PCLF.compute_polyhedral_pieces_pclf(
        f,
        graph,
        OPTIMIZER,
        Dict(v => PCLF.conic_partitions_2d(ORDER) for v in graph.verts);
        MLF = true,
    )
    isfinite(pclf.JSRapprox) || error("no certificate for $name")

    println(
        "\n=== ",
        name,
        ": certificate rate ",
        round(pclf.JSRapprox; digits = 6),
        " ===",
    )

    # The quotient is the expensive half -- minutes each -- and redrawing the figure must not pay for
    # it again. Delete the `.jld2` (or pass `force = true` below) after changing the graph, the
    # tolerances or `max_slices`, or a stale quotient is silently loaded against a fresh certificate.
    quotient, D = if isfile(cache_file)
        println("reusing the cached quotient from ", basename(cache_file))
        opt = import_optimizer_jld2(cache_file)
        MOI.get(opt, MOI.RawOptimizerAttribute("bisimulation_quotient")),
        MOI.get(opt, MOI.RawOptimizerAttribute("D"))
    else
        println("building the quotient ...")
        # `print_level = 1` on purpose: each build takes minutes, and the refinement counter is the
        # only progress signal there is. Silencing it makes a long run indistinguishable from a hung
        # one -- which cost a 13-minute run in this folder before the output was turned back on.
        r = build_quotient(
            problem,
            pclf;
            atol = 1e-4,
            level_tol = 1e-2,
            max_slices = 30,
            print_level = 1,
        )
        export_optimizer_jld2(r.optimizer, cache_file)
        r.quotient, r.D
    end
    # The two synthesis runs and the three volume integrals are the long pole once the quotient is
    # cached -- minutes each, against seconds for the quotient to load -- and redrawing the figure
    # must not pay for them either. They are small enough (three integer vectors and three scalars)
    # to sit in their own file, keyed to the quotient cache beside it.
    answers_file = replace(cache_file, ".jld2" => "_answers.jld2")
    answers = if isfile(answers_file)
        println("reusing the cached answers from ", basename(answers_file))
        jldopen(answers_file, "r") do file
            return Dict(k => file[k] for k in keys(file))
        end
    else
        println("solving their formula, both semantics ...")
        regions = Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3)
        ap_to_obs = Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3)

        exists = certified_set(f, quotient, regions, ap_to_obs).set
        # Both semantics are run for `∀`. `:language` is the answer this experiment wants; the
        # `:arbitrary` companion is what the code did before the distinction existed, and keeping it
        # side by side is what makes the difference a measurement rather than a claim.
        forall_env = ST.with_switching(f, HybridSystems.AutonomousSwitching())
        forall_lang = certified_set(forall_env, quotient, regions, ap_to_obs)
        forall_arb = certified_set(
            forall_env,
            quotient,
            regions,
            ap_to_obs;
            semantics = :arbitrary,
        )
        forall = forall_lang.set
        println(
            "  ∀ completions -- :language ",
            forall_lang.completions,
            " (genuine escapes)   :arbitrary ",
            forall_arb.completions,
            " (escapes + language-forbidden)",
        )

        # VOLUME, not cell count. The two quotients have different cells of different sizes covering
        # different regions, so neither raw counts nor fractions of cells compare across them -- that
        # is the cross-certificate warning in the README, and ignoring it made an earlier run of this
        # script report the two semantics moving in opposite directions purely as an artefact of cell
        # size. `get_volume` reports per-node volumes, so the total is their sum.
        vol(set) = sum(PCQ.get_volume(quotient, set; backend = CDDLib.Library()))

        a = Dict{String, Any}(
            "exists" => collect(exists),
            "forall" => collect(forall),
            "forall_arbitrary" => collect(forall_arb.set),
            "completions_language" => forall_lang.completions,
            "completions_arbitrary" => forall_arb.completions,
            "vol_exists" => vol(exists),
            "vol_forall" => vol(forall),
            "vol_forall_arbitrary" => vol(forall_arb.set),
            "vol_total" => vol(collect(keys(quotient.states))),
        )
        jldopen(answers_file, "w") do file
            for (k, v) in a
                file[k] = v
            end
        end
        a
    end

    return (;
        name,
        rate = pclf.JSRapprox,
        quotient,
        cells = length(quotient.states),
        exists = answers["exists"],
        forall = answers["forall"],
        forall_arbitrary = answers["forall_arbitrary"],
        completions_language = answers["completions_language"],
        completions_arbitrary = answers["completions_arbitrary"],
        vol_exists = answers["vol_exists"],
        vol_forall = answers["vol_forall"],
        vol_forall_arbitrary = answers["vol_forall_arbitrary"],
        vol_total = answers["vol_total"],
        slices = PCQ.num_slices(quotient),
        by_slice = PCQ.states_by_slice(quotient),
        degree = PCQ.outgoing_degree_stats(quotient),
    )
end

"""
One panel per (language, semantics): the certified region, drawn without axis limits.

This is the instrument the comparison actually needs. Cell counts are not comparable across two
quotients built from different certificates, but the regions they certify are -- they live in the
same state space, so the eye and the volume agree.

No limits are set here on purpose. A fixed ±10 window was wrong: the certificate's sublevel sets
reach far past the ±10 working set -- the covered region is a star, not a box -- so it showed only
the middle of each panel. `align_panels!` widens all four to a common window afterwards.
"""
function certified_panel(r, set, title)
    fig = plot(; aspect_ratio = :equal, legend = false, title = title, titlefontsize = 9)
    losing = collect(setdiff(Set(keys(r.quotient.states)), Set(set)))
    # Contours ON, which looks backwards and is not. The white haze over these panels is not strokes
    # -- the recipe already forces `linealpha = 0` when `show_contours` is false -- it is GR leaving
    # an anti-aliased seam between each pair of adjacent filled polygons, and at 10 611 cells the
    # seams outweigh the fill. Turning contours on strokes every cell in `linecolor := c`, its OWN
    # fill colour, so the strokes are invisible except that they cover the seams. `linewidth = 0`
    # cannot fix this; it makes it worse.
    _plot_winning_losing!(
        fig,
        r.quotient,
        collect(set),
        losing,
        nothing;
        show_contours = true,
        linewidth = 0.7,
    )
    plot!(
        fig,
        problem;
        plot_region = false,
        observation_region_alpha = 0.0,
        observation_colors = OBSERVATION_COLORS,
        observation_linewidth = 2.0,
    )
    return fig
end

"""
Put every panel on one window: the union of what Plots auto-scaled each of them to.

Reading the limits back is what makes this cheap. The alternative -- bounding the cells directly --
means a support-function LP per part over some 12 000 cells, where each panel already knows its own
extent from the shapes it just drew. The union matters because the two quotients cover differently
sized stars, and on separate windows two regions of different size look alike.
"""
function align_panels!(panels)
    xs = reduce(vcat, [collect(Plots.xlims(p)) for p in panels])
    ys = reduce(vcat, [collect(Plots.ylims(p)) for p in panels])
    lims = (extrema(xs), extrema(ys))
    for p in panels
        xlims!(p, lims[1]...)
        ylims!(p, lims[2]...)
    end
    return lims
end

"""
Cells per slice, and the growth factor between consecutive slices.

The quotient's size is not decided by the contraction rate alone. `thm:partition-growth` bounds the
cascade by the graph's out-degree -- p(i+1) <= p(i) + r * prod over outgoing edges of p(i) -- so a
better rate buys fewer slices while the out-degree decides how fast each slice grows. Both matter and
the profile shows them separately: the number of slices is the depth, the growth factor is the
branching. It is also the quantity that governs online cost, since a controller locating a state
walks this structure.
"""
function slice_profile(r)
    println("\n", r.name, ": ", r.slices, " slices, out-degree ", r.degree)
    counts = [get(Dict(r.by_slice), i, 0) for i in 1:(r.slices)]
    println(rpad("  slice", 10), rpad("cells", 9), "growth vs previous")
    for (i, c) in enumerate(counts)
        ratio =
            i == 1 || counts[i - 1] == 0 ? "-" :
            string(round(c / counts[i - 1]; digits = 2))
        println(rpad("  $i", 10), rpad(c, 9), ratio)
    end
    return nothing
end

# Unconstrained: the De Bruijn graph of order 1, which is E8's certificate.
unconstrained = both_semantics(
    "unconstrained",
    PCLF.generate_DeBruijn_edges(2, 1),
    joinpath(@__DIR__, "exp1_unconstrained.jld2"),
)

# Constrained: no two consecutive 1s.
no_11 = PCLF.edgeList_to_LabDigraph([(0, 1, 1), (0, 0, 2), (1, 0, 2)])
constrained = both_semantics("no 11", no_11, joinpath(@__DIR__, "exp1_no11.jld2"))

println("\nGol-Lazar-Belta, unmodified, their formula. Conic order ", ORDER, ".\n")
println(
    rpad("language", 16),
    rpad("rate", 11),
    rpad("cells", 8),
    rpad("∃ certified", 13),
    "∀ verified",
)
for r in (unconstrained, constrained)
    println(
        rpad(r.name, 16),
        rpad(round(r.rate; digits = 6), 11),
        rpad(r.cells, 8),
        rpad(length(r.exists), 13),
        length(r.forall),
    )
end

println(
    """

Prediction: constraining the switching should SHRINK the existential set (the controller loses
options) and GROW the universal one (the adversary loses options). Fractions are the honest
instrument here, since the two quotients have different cell counts.""",
)
for r in (unconstrained, constrained)
    slice_profile(r)
end

for r in (unconstrained, constrained)
    println(
        rpad(r.name, 16),
        "∃ ",
        rpad(round(100 * length(r.exists) / r.cells; digits = 1), 8),
        "%   ∀ ",
        round(100 * length(r.forall) / r.cells; digits = 1),
        "%",
    )
end

# ---------------------------------------------------------
# Volumes: the only cross-quotient comparison that is valid
# ---------------------------------------------------------

println("\nCertified VOLUME (cell counts are not comparable across certificates):\n")
println(
    rpad("language", 16),
    rpad("covered", 12),
    rpad("∃ volume", 12),
    rpad("∃ % of covered", 16),
    rpad("∀ volume", 12),
    "∀ % of covered",
)
for r in (unconstrained, constrained)
    println(
        rpad(r.name, 16),
        rpad(round(r.vol_total; digits = 2), 12),
        rpad(round(r.vol_exists; digits = 2), 12),
        rpad(round(100 * r.vol_exists / r.vol_total; digits = 1), 16),
        rpad(round(r.vol_forall; digits = 2), 12),
        round(100 * r.vol_forall / r.vol_total; digits = 1),
    )
end

# ---------------------------------------------------------
# The completion, which is what decided the old ∀ answer
# ---------------------------------------------------------

println("\nPessimistic completion, and what it was charging the adversary for:\n")
println(
    rpad("language", 16),
    rpad("nodes restrict?", 17),
    rpad(":language", 12),
    rpad(":arbitrary", 12),
    rpad("∀ (:language)", 15),
    "∀ (:arbitrary)",
)
for r in (unconstrained, constrained)
    println(
        rpad(r.name, 16),
        rpad(r.name == "no 11" ? "yes" : "no", 17),
        rpad(r.completions_language, 12),
        rpad(r.completions_arbitrary, 12),
        rpad(length(r.forall), 15),
        length(r.forall_arbitrary),
    )
end
println(
    """

On a complete graph the two semantics coincide and the columns agree. Under no-11 the difference is
every node-1 cell: `:arbitrary` credits the adversary with a mode-1 move the plant cannot make, and
because `∀` fails a state when ANY successor fails, that loss propagates backward and empties the
answer. `:language` completes only what genuinely escaped the covered region.""",
)

gr()
panels = [
    certified_panel(unconstrained, unconstrained.exists, "unconstrained, ∃ (synthesis)"),
    certified_panel(unconstrained, unconstrained.forall, "unconstrained, ∀ (verification)"),
    certified_panel(constrained, constrained.exists, "no 11, ∃ (synthesis)"),
    certified_panel(constrained, constrained.forall, "no 11, ∀ (verification)"),
]
lims = align_panels!(panels)
println(
    "\npanel window: x ∈ [",
    round(lims[1][1]; digits = 1),
    ", ",
    round(lims[1][2]; digits = 1),
    "], y ∈ [",
    round(lims[2][1]; digits = 1),
    ", ",
    round(lims[2][2]; digits = 1),
    "]",
)
fig = plot(panels...; layout = (2, 2), size = (900, 900))
savefig(fig, joinpath(@__DIR__, "fig_constraint_both_semantics.png"))
println("\nwrote fig_constraint_both_semantics.png")
