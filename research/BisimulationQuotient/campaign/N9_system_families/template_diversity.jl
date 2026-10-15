# Can an incomplete graph carry distinct per-node templates?
#
# Mechanism A (fewer cells) needs an INCOMPLETE graph; the facet advantage needs DISTINCT pieces.
# Every family in the folder has exactly one of the two: the rotation family has distinct pieces on
# a complete graph, the dual De Bruijn an incomplete graph with identical pieces. A family with both
# should win on cells and facets at once, which is the largest speedup available.
#
# First attempt, using the `:alternating` template mode -- and it fails:
#
#   θ        k=0 identical   k=1 dual identical   k=1 dual alt   k=2 dual alt
#   π/12     0.734833        0.734833             0.761890       0.794147
#   π/6      0.717560        0.717566             INFEASIBLE     INFEASIBLE
#   π/4      0.742542        0.742542             INFEASIBLE     INFEASIBLE
#   π/3      0.758582        0.758588             INFEASIBLE     INFEASIBLE
#
# But the test is badly designed and this is NOT a refutation. `:alternating` gives one node the
# rotated template and the other the IDENTITY, so it does not diversify the two nodes, it degrades
# one of them. And mismatched templates are strictly harder to satisfy along an edge, since
# V_d(A_m x) ≤ γ V_s(x) couples them -- an over-constrained conic program is what one should expect.
#
# The right test diversifies without degrading: give node s the template R(s·θ), distinct and
# equally expressive, which is precisely the structure that makes the ℤ_q family of N5 work. That
# needs a `:phase` mode in `rotation_templates`; until it is run, the combination is untested.
#
# Incidental confirmation: the dual De Bruijn's rate equals the single node's at every θ, so memory
# buys nothing in certificate quality on this system and the 174-vs-379 cell reduction of N1 really
# is purely structural.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

# The README's loosened tolerances: Clarabel's defaults do not converge reliably on the larger
# conic partitions, and report a certificate they cannot find as `Inf` rather than as a bound.
const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 4000,
    "tol_feas" => 1e-6,
    "tol_gap_abs" => 1e-6,
    "tol_gap_rel" => 1e-6,
)

function rate_of(f, k, dual, mode, θ)
    graph = PCLF.generate_DeBruijn_edges(2, k; dual = dual)
    nodes = sort(collect(graph.verts); by = string)
    return PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
        f,
        graph,
        OPTIMIZER;
        Gmats = rotation_templates(nodes; θ = θ, mode = mode),
        MLF = true,
        verbose = false,
    ).JSRapprox
end

function safe_rate(f, k, dual, mode, θ)
    value = try
        rate_of(f, k, dual, mode, θ)
    catch
        Inf
    end
    return isfinite(value) ? string(round(value; digits = 6)) : "INFEASIBLE"
end

(; f) = two_mode_problem()

configurations = [
    ("k=0 ident", 0, false, :rotation),
    ("k=1 dual ident", 1, true, :rotation),
    ("k=1 dual ALT", 1, true, :alternating),
    ("k=2 dual ALT", 2, true, :alternating),
]

println("Does an incomplete graph admit distinct per-node templates?")
println(rpad("θ", 10), join([rpad(name, 18) for (name, _, _, _) in configurations]))
for θ in [π / 12, π / 6, π / 4, π / 3]
    cells = [rpad(round(θ; digits = 3), 10)]
    for (_, k, dual, mode) in configurations
        push!(cells, rpad(safe_rate(f, k, dual, mode, θ), 18))
    end
    println(join(cells))
end
