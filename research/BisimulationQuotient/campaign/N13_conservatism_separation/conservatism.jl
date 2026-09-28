# A tighter certificate at equal per-node budget.
#
# Both certificates exist; one is strictly tighter. This is a CONSERVATISM separation and must not
# be confused with the FEASIBILITY separation of N2, where the single node certifies nothing at all.
# They are different quantities, and conflating them is an overclaim a referee will catch.
#
# Measured on the observer-graph problem, one 4-cone conic template per node either way:
#
#   JSR lower bound (periodic, L <= 10), word 211        0.868807
#   4-node observer graph                                0.868964  -> 0.02 % above the bound
#   1 node, best of 41 orientations                      0.873586  -> 0.55 % above the bound
#   gap                                                  0.004622
#
# The identity orientation alone gives 0.9137. Quoting that was how this separation was first
# overstated by a factor of ten -- the third time in this work an under-searched baseline has done
# so, which is why the best-so-far curve is printed.
#
# So the path-complete family is essentially tight against the JSR while the single node is not.
# The four pieces are genuinely distinct (pairwise support-function distance 0.025 to 0.470), the
# graph is incomplete, and the s.m.p. is the length-3 word 211 -- there is word structure to
# remember, which is exactly what N3 predicts.
#
# The single-node side is searched over many template orientations, with best-so-far reported
# against the search budget. Under-searching a baseline has manufactured separations twice in this
# work; it must not happen a third time.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Random
using LinearAlgebra

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 3000,
    "verbose" => false,
    "tol_feas" => 1e-6,
    "tol_gap_abs" => 1e-6,
    "tol_gap_rel" => 1e-6,
)

"""
Best rate a single node certifies at the shared template, as a function of the search budget.

Returns the running best after each draw, so that under-searching is visible rather than hidden in
a single number. The identity orientation is always tried first.
"""
function single_node_curve(f, partition; draws = 40, seed = 11)
    graph = PCLF.edgeList_to_LabDigraph([(1, 1, 1), (1, 1, 2)])
    rng = Random.MersenneTwister(seed)
    best = Inf
    curve = Float64[]
    for i in 0:draws
        M = i == 0 ? Matrix{Float64}(I, 2, 2) : Matrix(qr(randn(rng, 2, 2)).Q)
        rate = try
            PCLF.compute_polyhedral_pieces_pclf(
                f,
                graph,
                OPTIMIZER,
                Dict(v => [M * c for c in partition] for v in graph.verts);
                MLF = true,
            ).JSRapprox
        catch
            Inf
        end
        best = min(best, rate)
        push!(curve, best)
    end
    return curve
end

spectral_radius(M) = maximum(abs.(LinearAlgebra.eigvals(M)))

function jsr_lower_bound(A, Lmax)
    best, best_word = 0.0, Int[]
    for L in 1:Lmax, code in 0:(2 ^ L - 1)
        word = [((code >> (i - 1)) & 1) + 1 for i in 1:L]
        P = Matrix{Float64}(LinearAlgebra.I, 2, 2)
        for m in word
            P = A[m] * P
        end
        r = spectral_radius(P)^(1 / L)
        if r > best
            best, best_word = r, word
        end
    end
    return best, best_word
end

(; f) = observer_graph_problem(; p = 1.7)
partition = PCLF.base_conic_partition_2d()

bound, word = jsr_lower_bound(UT.mode_matrices(f), 10)
pclf_rate = observer_graph_pclf(f).JSRapprox
curve = single_node_curve(f, partition)

println("Conservatism separation on the observer-graph problem, 4 cones per node")
println(
    "  JSR lower bound (periodic, L ≤ 10) : ",
    round(bound; digits = 6),
    "  at word ",
    word,
)
println(
    "  4-node observer graph              : ",
    round(pclf_rate; digits = 6),
    "   (",
    round(100 * (pclf_rate / bound - 1); digits = 2),
    " % above the bound)",
)
println(
    "  1 node, best of ",
    length(curve),
    " orientations   : ",
    round(curve[end]; digits = 6),
    "   (",
    round(100 * (curve[end] / bound - 1); digits = 2),
    " % above the bound)",
)
println(
    "  gap                                : ",
    round(curve[end] - pclf_rate; digits = 6),
)

println(
    "\n  single-node best-so-far against search budget (under-searching must be visible):",
)
for i in [1, 2, 5, 10, 20, length(curve)]
    i <= length(curve) &&
        println("     after ", rpad(i, 4), " draws : ", round(curve[i]; digits = 6))
end
