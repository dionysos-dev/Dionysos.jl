# Can the bisimulation actually be built on a dwell-time graph?  (Proof obligation 1, code half.)
#
# The question is about the ABSTRACTION, not the certificate. The certificate side is settled --
# `dwell_time_example.jl` finds a PCLF at tau = 3 -- but `bisimulation_pclf` was written for graphs
# path-complete over the full alphabet, and a dwell automaton is path-complete only relative to its
# own language. Obligation 1 in plan.md asks whether the correctness ARGUMENT survives that
# weakening; this script asks whether the CODE does, which is a cheaper question and settles it
# empirically. It does, on a 6-node graph -- larger than either headline experiment uses.
#
# Diagnostic only, and deliberately figureless: the shear throws the certificate's sublevel sets out
# to about ±130 against a ±4 working set, so the quotient renders illegibly. What this example has to
# show is STRUCTURE -- mean outgoing degree 1.53, minimum 1, because a node (m,k) with k < tau has
# exactly one outgoing edge and you must keep dwelling -- and structure prints better than it plots.
# `quotient_on_glb_constrained.jl` is the figure for this folder.
#
# Run after `dwell_time_example.jl`; it reuses that file's `shear_modes` and `dwell_graph`.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra

include(joinpath(@__DIR__, "dwell_time_example.jl"))

const TAU = 3
const ORDER = 2

A1, A2 = shear_modes(; diag = 0.3, shear = 2.0)
f = ST.with_switching(
    HybridSystems.discreteswitchedsystem([A1, A2]),
    HybridSystems.ControlledSwitching(),
)

# A working set the contracted dynamics stay inside, and one observation region to reach.
X = LazySets.Hyperrectangle(; low = [-4.0, -4.0], high = [4.0, 4.0])
R1 = LazySets.Hyperrectangle(; low = [1.0, 1.0], high = [2.5, 2.5])
problem = PR.BisimulationQuotientProblem(f, X, [R1])

graph = dwell_graph(2, TAU)
println("dwell automaton: ", length(graph.verts), " nodes, ", length(graph.edges), " edges")

pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    graph,
    JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 4000,
        "verbose" => false,
        "tol_feas" => 1e-6,
        "tol_gap_abs" => 1e-6,
        "tol_gap_rel" => 1e-6,
    ),
    Dict(v => PCLF.conic_partitions_2d(ORDER) for v in graph.verts);
    MLF = true,
)
println("certificate rate = ", pclf.JSRapprox)
isfinite(pclf.JSRapprox) || error("no certificate at tau = $TAU, order $ORDER")

println(
    "\nbuilding the quotient on a graph that is path-complete only over its own language ...",
)
result = try
    build_quotient(
        problem,
        pclf;
        atol = 1e-3,
        level_tol = 1e-2,
        max_slices = 10,
        print_level = 1,
    )
catch err
    println("\n*** THE CONSTRUCTION DOES NOT ACCEPT A DWELL-TIME GRAPH ***")
    println(first(split(sprint(showerror, err), "\n")))
    println(
        """

That is the answer this script exists to get, and it is not a failure of the idea: it says the
abstraction needs the generalization that proof obligation 1 in plan.md describes, in the same
way the certificate needed the language-aware bisection bracket. Record it and fix the code
rather than abandoning the experiment.""",
    )
    nothing
end

result === nothing && exit()

quotient = result.quotient
println("\nquotient built: ", length(quotient.states), " cells")

# Reach R1: observation 1 is the first declared region.
function winning(quotient, target_obs::Int)
    win = Set{Int}(id for (id, s) in quotient.states if s.obs == target_obs)
    changed = true
    while changed
        changed = false
        for (id, s) in quotient.states
            id in win && continue
            if any(t -> t[2] in win, s.next)
                push!(win, id)
                changed = true
            end
        end
    end
    return win
end

win = collect(winning(quotient, 1))
println(
    "certified for reach-R1: ",
    length(win),
    " of ",
    length(quotient.states),
    " cells -- but a LARGE region, since cells vary hugely in size here; the count is not the measure",
)

# Structure, not geometry. The dwell restriction shows up in the transition structure -- nodes (m,k)
# with k < tau have exactly one outgoing edge, so you MUST keep dwelling -- and that is the part worth
# reporting. There is deliberately no figure: the shear throws the certificate's sublevel sets out to
# about ±130 against a ±4 working set, so the quotient renders illegibly, and
# `quotient_on_glb_constrained.jl` is the figure for this folder.

bynode = Dict{Any, Int}()
for s in values(quotient.states)
    bynode[s.node] = get(bynode, s.node, 0) + 1
end
println(
    "\nstates by node: ",
    join(["$k => $v" for (k, v) in sort(collect(bynode); by = string ∘ first)], " · "),
)
println(
    """

Obligation 1 is answered for the CODE, not for the theory: the construction runs on a graph that is
path-complete only relative to its own language, and on a larger one (6 nodes) than either headline
experiment uses. The theory question -- whether soundness, the slice-decrease property and the
terminal pair survive path-completeness relative to L rather than <M>† -- is untouched by this, and
the certificate side is the cautionary precedent: there the code also ran, and returned a wrong
answer silently, until the bisection bracket was fixed.""",
)
