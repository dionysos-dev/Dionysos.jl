# Does the spectrum-maximizing product predict whether memory can help?
#
# Memory encodes which words have been seen. If the extremal behaviour of the system is a single
# mode repeated, there is nothing to remember, and a path-complete family cannot beat a common one
# at any template budget. The length of the s.m.p. is therefore a cheap a-priori test, costing a
# handful of eigenvalue computations, of whether the method can help on a given system.
#
# Measured result:
#
#   system            JSR lb    attained by     max_m ρ(A_m)   ratio      memory gap
#   gol_lazar_belta   0.85585   2111111         0.85580        1.00006    0.0
#   observer_graph    0.86881   211             0.80016        1.086      0.0046
#   two_mode          0.70395   2111111111      0.70000        1.00564    0.0
#
# THE PREDICTOR IS THE RATIO, NOT THE WORD LENGTH. `two_mode` has a length-10 s.m.p. word and yet
# memory buys nothing on it -- N1 measures the dual De Bruijn and the single node certifying the
# identical rate 0.717572 -- because the word beats the best single mode by only 0.56 %. A long word
# whose product barely exceeds a single mode carries no structure worth remembering.
#
# So the rule is: when max_w ρ(A_w)^{1/|w|} is close to max_m ρ(A_m), the extremal behaviour is
# essentially "stay in one mode", and no path-complete family can beat a common one. Three systems,
# three consistent calls. On Gol-Lazar-Belta this is confirmed independently by the budget sweep
# below, where one, two and four nodes certify at identical rates at conic orders 1, 2 and 3 alike.
#
# The same reasoning explains an earlier negative result: for random
# planar systems a four-facet polytope is essentially an extremal norm, so the s.m.p. is short and
# the common certificate is already tight.
#
# Caveat to state wherever this is used: a periodic lower bound over words of length ≤ L is a lower
# bound, not the JSR. An s.m.p. longer than L is invisible to it, so report L and the best word.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra

spectral_radius(M) = maximum(abs.(LinearAlgebra.eigvals(M)))

"""
Best periodic rate `ρ(A_w)^{1/|w|}` over words of length at most `Lmax`, with the word attaining it.

A lower bound on the joint spectral radius, and the instrument this file is built on.
"""
function smp_lower_bound(A, Lmax::Int)
    best, best_word = 0.0, Int[]
    M = length(A)
    for L in 1:Lmax, code in 0:(M ^ L - 1)
        word = [(div(code, M^(i - 1)) % M) + 1 for i in 1:L]
        P = Matrix{Float64}(LinearAlgebra.I, size(A[1])...)
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

function predictor_row(name, f; Lmax = 10)
    A = UT.mode_matrices(f)
    bound, word = smp_lower_bound(A, Lmax)
    best_single = maximum(spectral_radius, A)
    println(
        rpad(name, 18),
        rpad(round(bound; digits = 5), 10),
        rpad(join(word), 14),
        rpad(round(best_single; digits = 5), 12),
        round(bound / best_single; digits = 5),
    )
    return nothing
end

"""
Rate certified at each conic order by the De Bruijn graph of order `k`, on one system.

The budget sweep: if memory buys nothing at *every* template budget, the system is structurally
memory-free and no benchmark built on it can demonstrate the mechanism.
"""
function budget_sweep(name, f; orders = 1:3, ks = 0:2)
    optimizer = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 2000,
        "tol_feas" => 1e-6,
        "tol_gap_abs" => 1e-6,
        "tol_gap_rel" => 1e-6,
    )
    println("\n", name, ": rate by (template budget, graph order)")
    println(
        rpad("conic order", 13),
        rpad("rows/node", 11),
        join([rpad("k=$k ($(2^k) node$(2^k == 1 ? "" : "s"))", 22) for k in ks]),
    )
    for order in orders
        partition = PCLF.conic_partitions_2d(order)
        cells = [rpad(string(order), 13), rpad(string(length(partition)), 11)]
        for k in ks
            graph = PCLF.generate_DeBruijn_edges(2, k)
            rate = try
                string(
                    round(
                        PCLF.compute_polyhedral_pieces_pclf(
                            f,
                            graph,
                            optimizer,
                            Dict(v => partition for v in graph.verts);
                            MLF = true,
                        ).JSRapprox;
                        digits = 6,
                    ),
                )
            catch
                "FAILED"
            end
            push!(cells, rpad(rate, 22))
        end
        println(join(cells))
    end
    return nothing
end

systems = [
    ("gol_lazar_belta", gol_lazar_belta_problem().f),
    ("observer_graph", observer_graph_problem(; p = 1.7).f),
    ("two_mode", two_mode_problem().f),
]

println("The s.m.p. predictor: is the extremal behaviour a long word or a single mode?")
println(
    rpad("system", 18),
    rpad("JSR lb", 10),
    rpad("attained by", 14),
    rpad("max ρ(A_m)", 12),
    "ratio",
)
for (name, f) in systems
    predictor_row(name, f)
end

# The negative control: memory buys nothing on Gol-Lazar-Belta at any budget.
budget_sweep("gol_lazar_belta", gol_lazar_belta_problem().f)
