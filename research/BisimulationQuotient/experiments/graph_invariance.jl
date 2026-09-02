# The graph moves cost, not the answer.
#
# One system, one co-safe LTL specification, four path-complete graphs. The certified winning set
# is compared pointwise on a grid, alongside every cost measure. The certified set is invariant
# while the cost is not, so the graph is a design parameter that carries no conservatism penalty.
#
# The last graph is deliberately path-complete but *not* complete. For it the quotient is only
# behaviourally equivalent to the concrete system, so the usual bisimulation argument does not
# apply and only the synthesis-invariance result predicts agreement — making it the sharpest
# available test. It is also, on this instance, the cheapest of the four.
#
# Volumes are deliberately not used as the instrument: `max_slices` truncates the sublevel family
# at a different outer level for each certificate, so the graphs do not cover the same region and
# raw volumes differ for reasons unrelated to the graph. A pointwise test has no such confound.

include(joinpath(@__DIR__, "common.jl"))

using Spot
using Printf

const SPEC = ltl"F(R1 & F(D))"
const X0 = SVector(2.3, 1.5)
const ATOL = 1e-3
const MAX_SLICES = 8

# A grid over the state set, plus the initial state the other scripts use.
const PROBES = vcat([SVector(x, y) for x in -1.9:0.4:1.9 for y in -1.9:0.4:1.9], [X0])

"""
Whether `x` lies in some winning cell of `quotient`.

This is the existential semantics: a concrete state is controllable when *some* lifted copy of it
is. Under a complete graph the choice is immaterial, since all copies of a state are bisimilar.
"""
function is_winning(quotient, controllable_set, x)
    return any(qid -> x ∈ quotient.states[qid].set, controllable_set)
end

function run_graph(label, graph)
    (; f, problem, R1, R2) = two_mode_problem()
    nodes = sort(collect(graph.verts))

    pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
        f,
        graph,
        JuMP.optimizer_with_attributes(
            Clarabel.Optimizer,
            "max_iter" => 2000,
            "verbose" => false,
        );
        Gmats = rotation_templates(nodes),
        MLF = true,
        verbose = false,
    )
    isfinite(pclf.JSRapprox) || return nothing

    t0 = time()
    (; quotient, D) = build_quotient(
        problem,
        pclf;
        atol = ATOL,
        max_levels = 100,
        max_slices = MAX_SLICES,
        print_level = 0,
    )
    build_s = time() - t0

    t0 = time()
    result = synthesize_cosafe_ltl(
        f,
        quotient,
        Dionysos.spot_stepper(SPEC),
        Dict(:D => D, :R1 => R1, :R2 => R2),
        Dict(:D => -1, :R1 => 1, :R2 => 2),
        X0;
        print_level = 0,
    )
    synth_s = time() - t0

    _, faces = PCQ.cell_complexities(quotient)
    winning = [is_winning(quotient, result.controllable_set, x) for x in PROBES]

    return (;
        label,
        nodes = length(nodes),
        complete = PCLF.is_complete(graph, 1:2),
        cells = length(quotient.states),
        facets = sum(faces),
        max_facets = maximum(faces),
        build_s,
        synth_s,
        winning,
    )
end

const GRAPHS = [
    ("De Bruijn k=0", PCLF.generate_DeBruijn_edges(2, 0)),
    ("De Bruijn k=1", PCLF.generate_DeBruijn_edges(2, 1)),
    ("De Bruijn k=2", PCLF.generate_DeBruijn_edges(2, 2)),
    ("dual De Bruijn k=1", PCLF.generate_DeBruijn_edges(2, 1; dual = true)),
]

results = filter(!isnothing, [run_graph(label, graph) for (label, graph) in GRAPHS])

println()
println(
    rpad("graph", 21),
    rpad("nodes", 7),
    rpad("complete", 10),
    rpad("cells", 8),
    rpad("Σ facets", 10),
    rpad("max fct", 9),
    rpad("build s", 9),
    rpad("synth s", 9),
    "winning",
)
println("-"^96)
for r in results
    @printf(
        "%-21s%-7d%-10s%-8d%-10d%-9d%-9.1f%-9.2f%d/%d\n",
        r.label,
        r.nodes,
        r.complete ? "yes" : "no",
        r.cells,
        r.facets,
        r.max_facets,
        r.build_s,
        r.synth_s,
        count(r.winning),
        length(PROBES),
    )
end

reference = results[1].winning
println()
println("Pointwise agreement with $(results[1].label) over $(length(PROBES)) probes:")
for r in results
    disagreements = count(r.winning .!= reference)
    println(
        "  ",
        rpad(r.label, 21),
        disagreements == 0 ? "identical" : "$disagreements DISAGREEMENTS",
    )
end

cells = [r.cells for r in results]
facets = [r.facets for r in results]
println()
@printf(
    "Cost spread at a constant answer: cells %d–%d (%.1f×), facets %d–%d (%.1f×)\n",
    minimum(cells),
    maximum(cells),
    maximum(cells) / minimum(cells),
    minimum(facets),
    maximum(facets),
    maximum(facets) / minimum(facets),
)
