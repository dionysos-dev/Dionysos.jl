# E3: the running example when the plant never applies mode 1 twice in a row (paper, Section 10.3).
#
#     julia --project=. experiments/E3_constrained_plant.jl   [CACHE=0] [PROBES=n]
#
# The constraint belongs to the plant. The certificate on the constraint graph represents it, so its
# quotient answers synthesis and verification exactly. The quotient of the conference running
# example, built for arbitrary switching, covers the constraint without respecting it: its
# verification answer is sound but conservative and its synthesis answer is unsound. Both are
# evaluated at the same probe points, and the inclusions the theory predicts are checked there.
include(joinpath(dirname(@__DIR__), "common.jl"))

const PROBES = parse(Int, get(ENV, "PROBES", "211"))

f = switched([A1, A2]; controlled = true)
f_forall = switched([A1, A2]; controlled = false)
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])

# Constraint state 0: mode 1 may be applied. Constraint state 1: mode 1 was just applied.
constraint = PCLF.edgeList_to_LabDigraph([(0, 0, 2), (0, 1, 1), (1, 0, 2)])
debruijn = PCLFProducer.generate_DeBruijn_edges(2, 1)
init_constrained = Dict(0 => [0], 1 => [1])
init_free = Dict(c => sort(collect(debruijn.verts)) for c in (0, 1))

@assert PCLF.respects_constraint(constraint, init_constrained, constraint)
@assert PCLF.covers_constraint(constraint, init_constrained, constraint)
@assert PCLF.covers_constraint(debruijn, init_free, constraint)
@assert !PCLF.respects_constraint(debruijn, init_free, constraint)

println("="^84, "\nE3 — the running example under the constraint of Fig. 1(a)\n", "="^84)
free = certificate(f, debruijn; with_common = true)
cons = certificate(f, constraint)
Printf.@printf(
    "rate             arbitrary switching %.6f, constrained %.6f\n",
    free.γ,
    cons.γ
)
Printf.@printf("slices           arbitrary switching %d, constrained %d\n", free.N, cons.N)

built_free = build("E3_free", problem, free, free.pclf)
built_cons = build("E3_constrained", problem, cons, cons.pclf)
stats_free, stats_cons = cell_stats(built_free.quotient), cell_stats(built_cons.quotient)
Printf.@printf(
    "cells            arbitrary switching %d, constrained %d\n",
    stats_free.cells,
    stats_cons.cells
)
stats_free.cells == 4625 ||
    @warn "the arbitrary-switching quotient differs from the conference running example (4625 cells)" stats_free.cells

won = (;
    free_exists = solve(f, built_free),
    free_forall = solve(f_forall, built_free),
    cons_exists = solve(f, built_cons),
    cons_forall = solve(f_forall, built_cons),
)
# Timed once both quotients exist, as in the conference experiments.
seconds_free = build_seconds(problem, free, free.pclf)
seconds_cons = build_seconds(problem, cons, cons.pclf)
Printf.@printf(
    "build time       arbitrary switching %.1f s, constrained %.1f s\n",
    seconds_free,
    seconds_cons
)

points, area = probe_grid(PROBES)
cells_free = Dict(nd => locate(built_free.quotient, nd, points) for nd in debruijn.verts)
cells_cons = Dict(nd => locate(built_cons.quotient, nd, points) for nd in constraint.verts)

numbers = Pair{String, Any}[
    "EThreeRateFree" => fmt(free.γ; digits = 3),
    "EThreeRateCons" => fmt(cons.γ; digits = 3),
    "EThreeCellsFree" => fmt(stats_free.cells),
    "EThreeCellsCons" => fmt(stats_cons.cells),
    "EThreeN" => fmt(cons.N - 1),
    "EThreeTimeFree" => fmt_seconds(seconds_free),
    "EThreeTimeCons" => fmt_seconds(seconds_cons),
]
panels = Plots.Plot[]
for (c, tag) in ((0, "Zero"), (1, "One"))
    nodes_cons, nodes_free = init_constrained[c], init_free[c]
    domain = covered(cells_free, nodes_free) .& covered(cells_cons, nodes_cons)
    exists_exact =
        project(cells_cons, won.cons_exists, nodes_cons; quantifier = :exists) .& domain
    forall_exact =
        project(cells_cons, won.cons_forall, nodes_cons; quantifier = :forall) .& domain
    exists_free =
        project(cells_free, won.free_exists, nodes_free; quantifier = :exists) .& domain
    forall_free =
        project(cells_free, won.free_forall, nodes_free; quantifier = :forall) .& domain
    # Theorem 5.5: the exact sets are those of the constrained plant, the arbitrary-switching ones
    # those of the unconstrained plant, and fewer signals shrink ∃ and enlarge ∀.
    broken_exists = count(exists_exact .& .!exists_free)
    broken_forall = count(forall_free .& .!forall_exact)
    Printf.@printf(
        "constraint state %d  domain %.1f   ∃: exact %.1f, arbitrary switching %.1f   ∀: exact %.1f, arbitrary switching %.1f\n",
        c,
        area * count(domain),
        area * count(exists_exact),
        area * count(exists_free),
        area * count(forall_exact),
        area * count(forall_free),
    )
    Printf.@printf(
        "                    points violating the predicted inclusions: ∃ %d, ∀ %d\n",
        broken_exists,
        broken_forall
    )
    append!(
        numbers,
        [
            "EThreeDomain$tag" => fmt(area * count(domain)),
            "EThreeExistsExact$tag" => fmt(area * count(exists_exact)),
            "EThreeExistsFree$tag" => fmt(area * count(exists_free)),
            "EThreeForallExact$tag" => fmt(area * count(forall_exact)),
            "EThreeForallFree$tag" => fmt(area * count(forall_free)),
            "EThreeBroken$tag" => fmt(broken_exists + broken_forall),
        ],
    )
    # Where ignoring the constraint errs: the arbitrary-switching ∃-set claims points no admissible
    # signal wins, and its ∀-set misses points every admissible signal wins.
    exists_category = [
        d ? (e ? 3 : (u ? 2 : 1)) : 0 for
        (d, e, u) in zip(domain, exists_exact, exists_free)
    ]
    forall_category = [
        d ? (a ? 3 : (e ? 2 : 1)) : 0 for
        (d, a, e) in zip(domain, forall_free, forall_exact)
    ]
    # No titles: with an equal aspect ratio GR clips those of the top row; the caption names the panels.
    push!(
        panels,
        category_panel(PROBES, points, exists_category, [LOST, UNSOUND, WON]; title = ""),
    )
    push!(
        panels,
        category_panel(PROBES, points, forall_category, [LOST, MISSED, WON]; title = ""),
    )
end
write_numbers("E3_constrained_plant", numbers)
println("numbers written to results/numbers_E3_constrained_plant.tex")

figure = Plots.plot(panels...; layout = (2, 2), size = (860, 860), dpi = 200)
Plots.savefig(figure, joinpath(RESULTS, "fig_E3_constrained_plant.png"))
println("figure written to results/fig_E3_constrained_plant.png")
