# Does a polyhedral PCLF actually exist on the constraint automaton?
#
# This is the one real risk in N2. The periodic bound says the constrained system is stabilizable at
# rate >= 0.912 for c = 1.25, but that is a statement about the system, not about whether a template
# of a given order reaches it. If no certificate is found at conic order 2 or 3, the experiment needs
# a different c or a richer template, and it is much cheaper to learn that here than after building
# a quotient.
#
# The system: Sigma_c = {c*A1, A2} with Gol-Lazar-Belta's matrices. For c >= 1.17 the constant word
# 1^omega diverges, so no common Lyapunov function exists at any complexity and the predecessor's
# construction has no input at all -- established by one eigenvalue, not by a search.
#
# The constraint automaton forbids two consecutive uses of mode 1:
#
#     q0 --1--> q1        q0 --2--> q0        q1 --2--> q0
#
# and q1 has no mode-1 edge, which is exactly the prohibition. Note this graph is NOT path-complete
# over {1,2}* -- it is path-complete relative to the constrained language, which is the point, and
# the reason the correctness argument has to be re-read before this becomes an experiment.
#
# For contrast the single-node graph (both modes as self-loops) is also solved: it represents the
# unconstrained system, and it must FAIL for every c in the window. That failure is the separation.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra

const A1 = [-0.65 0.32; -0.42 -0.92]
const A2 = [0.65 0.32; -0.42 -0.92]

spectral_radius(M) = maximum(abs.(LinearAlgebra.eigvals(M)))

scaled_system(c) = ST.with_switching(
    HybridSystems.discreteswitchedsystem([c .* A1, Matrix(A2)]),
    HybridSystems.ControlledSwitching(),
)

"""
The two-state automaton generating the language with no two consecutive 1s.
"""
no_11_graph() = PCLF.edgeList_to_LabDigraph([(0, 1, 1), (0, 0, 2), (1, 0, 2)])

"""
The single node with both modes: the unconstrained system, which must fail in the window.
"""
unconstrained_graph() = PCLF.edgeList_to_LabDigraph([(1, 1, 1), (1, 1, 2)])

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 4000,
    "verbose" => false,
    "tol_feas" => 1e-6,
    "tol_gap_abs" => 1e-6,
    "tol_gap_rel" => 1e-6,
)

function rate(f, graph, order)
    partition = PCLF.conic_partitions_2d(order)
    value = try
        PCLF.compute_polyhedral_pieces_pclf(
            f,
            graph,
            OPTIMIZER,
            Dict(v => partition for v in graph.verts);
            MLF = true,
        ).JSRapprox
    catch
        Inf
    end
    return value
end

fmt(v) = isfinite(v) ? string(round(v; digits = 6)) : "none"

println("Sigma_c = {c*A1, A2}, Gol-Lazar-Belta's matrices.")
println(
    "A certificate on the no-11 automaton is what N2 needs; the single node must fail.\n",
)
println(
    rpad("c", 7),
    rpad("rho(cA1)", 11),
    rpad("order", 7),
    rpad("1 node (uncon.)", 17),
    rpad("no-11 automaton", 17),
    "verdict",
)
for c in [1.25, 1.30], order in 1:3
    f = scaled_system(c)
    one = rate(f, unconstrained_graph(), order)
    con = rate(f, no_11_graph(), order)
    verdict = if isfinite(con) && con < 1.0 && !(isfinite(one) && one < 1.0)
        "*** separation: constrained certifies, unconstrained cannot ***"
    elseif isfinite(con) && con < 1.0
        "both certify -- c is too small, no separation"
    else
        "no certificate on the constrained language at this order"
    end
    println(
        rpad(c, 7),
        rpad(round(spectral_radius(c * A1); digits = 4), 11),
        rpad(order, 7),
        rpad(fmt(one), 17),
        rpad(fmt(con), 17),
        verdict,
    )
end
