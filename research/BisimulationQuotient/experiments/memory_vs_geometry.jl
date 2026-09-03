# Memory buys an abstraction that bounded geometry cannot.
#
# A family for which the single-node (common Lyapunov) construction admits no certificate at a
# given per-node facet budget, while a path-complete family at the *same* budget does — so the
# predecessor method returns no abstraction at all and this one returns a quotient.
#
#     A₁ = ρ R(2π/q),   A₂ = ρ R(-2π/q)
#     graph  ℤ_q  with edges (s, s+1, mode 1) and (s, s-1, mode 2)   -- complete, deterministic
#     PCLF   V_s(x) = ‖R^{-s} x‖_∞                                   -- four facets per node, any q
#
# Why it separates. Along the edge (s, s+1, 1),
#     V_{s+1}(A₁x) = ‖R^{-(s+1)} ρ R x‖_∞ = ρ ‖R^{-s}x‖_∞ = ρ V_s(x),
# and symmetrically for mode 2, so the cycle certifies at rate exactly ρ for every q with four
# facets per node. A *common* function with unit ball B must satisfy R^{±1}B ⊆ (γ/ρ)B; as γ → ρ
# this forces B invariant under the order-q rotation group, and a polytopic B then needs at least
# q facets. Restricted to four, the achievable rate is bounded away from ρ: for q = 5 the gauge of
# R_{72°}(1,1) is 1.260, so the best rate is ≈ 1.26ρ and any ρ > 1/1.26 ≈ 0.794 defeats it.
#
# This instantiates the facet-budget bound: bounding the per-node geometry at ℓ forces |S| to grow,
# because the total budget is bounded below by a quantity depending only on the system.
#
# The certificate is written down analytically rather than solved for, so no optimizer enters the
# claim; `verify_pclf_rate` checks it by sampling. The single-node side is searched over many
# template orientations, since under-searching it would manufacture the separation.

include(joinpath(dirname(@__DIR__), "common.jl"))

using Random
using Printf

rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

"""
Two modes rotating by ±2π/q and contracting by `ρ`.
"""
rotation_system(ρ, q) =
    HybridSystems.discreteswitchedsystem([ρ .* rotation(2π / q), ρ .* rotation(-2π / q)])

"""
The cycle graph on `ℤ_q` and its analytic PCLF: node `s` carries the square `‖R^{-s}x‖_∞ ≤ γ`.
"""
function cycle_pclf(ρ, q)
    θ = 2π / q
    edges = vcat(
        [(s, mod(s + 1, q), 1) for s in 0:(q - 1)],
        [(s, mod(s - 1, q), 2) for s in 0:(q - 1)],
    )
    graph = PCLF.edgeList_to_LabDigraph(edges)
    pieces = Dict{Any, PCLF.AbstractPiece}(
        s => PCLF.PolyhedralPiece(Matrix(rotation(-s * θ)), ones(2)) for s in 0:(q - 1)
    )
    return graph, PCLF.PCLF(graph, pieces, ρ)
end

"""
Largest observed `V_d(A_m x) / V_s(x)` over the graph edges, sampled at random points.

An independent check on the analytic certificate: it should return `ρ`, and it involves no solver.
"""
function verify_pclf_rate(ρ, q; nsample = 4000, seed = 3)
    θ = 2π / q
    modes = [ρ .* rotation(2π / q), ρ .* rotation(-2π / q)]
    V(s, x) = maximum(abs.(rotation(-s * θ) * x))
    rng = Random.MersenneTwister(seed)
    worst = 0.0
    for _ in 1:nsample
        x = randn(rng, 2)
        for s in 0:(q - 1)
            vs = V(s, x)
            vs < 1e-12 && continue
            worst = max(worst, V(mod(s + 1, q), modes[1] * x) / vs)
            worst = max(worst, V(mod(s - 1, q), modes[2] * x) / vs)
        end
    end
    return worst
end

"""
Best rate a single node can certify with a four-facet template, over `draws` orientations.

The identity is always included: a purely random search can miss the axis-aligned template and
would then understate what one node achieves.
"""
function single_node_rate(ρ, q; draws = 40, seed = 11)
    f = rotation_system(ρ, q)
    graph = PCLF.generate_DeBruijn_edges(2, 0)
    nodes = sort(collect(graph.verts))
    optimizer = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 3000,
        "verbose" => false,
        "tol_feas" => 1e-6,
        "tol_gap_abs" => 1e-6,
        "tol_gap_rel" => 1e-6,
    )
    rng = Random.MersenneTwister(seed)
    best = Inf
    for i in 0:draws
        M = i == 0 ? Matrix{Float64}(LinearAlgebra.I, 2, 2) : Matrix(qr(randn(rng, 2, 2)).Q)
        rate = try
            PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
                f,
                graph,
                optimizer;
                Gmats = Dict(v => M for v in nodes),
                MLF = true,
                verbose = false,
            ).JSRapprox
        catch
            Inf
        end
        best = min(best, rate)
    end
    return best
end

# ---------------------------------------------------------
# The separation, over a range of (q, ρ)
# ---------------------------------------------------------

println(
    "A₁ = ρR(2π/q), A₂ = ρR(-2π/q); cycle graph on ℤ_q. Four facets per node throughout.",
)
println()
println(
    rpad("q", 4),
    rpad("ρ", 7),
    rpad("one node", 12),
    rpad("cycle", 10),
    rpad("verified", 11),
    "verdict",
)
println("-"^62)

for q in (5, 7), ρ in (0.80, 0.88, 0.94)
    _, pclf = cycle_pclf(ρ, q)
    verified = verify_pclf_rate(ρ, q)
    one_node = single_node_rate(ρ, q)
    separates = one_node >= 1 && verified < 1
    @printf(
        "%-4d%-7.2f%-12s%-10.5f%-11.5f%s\n",
        q,
        ρ,
        isfinite(one_node) ? @sprintf("%.5f", one_node) : "Inf",
        pclf.JSRapprox,
        verified,
        separates ? "SEPARATION" :
        (one_node < 1 ? "one node suffices" : "neither certifies"),
    )
end

# ---------------------------------------------------------
# From certificate to abstraction
# ---------------------------------------------------------
# The certificate gap is only the premise. What matters is that one construction yields a
# bisimulation quotient and the other cannot be started.

const Q, RHO = 5, 0.88
graph, pclf = cycle_pclf(RHO, Q)

println()
println("Building the quotient for q = $Q, ρ = $RHO")
println(
    "  graph: complete = ",
    PCLF.is_complete(graph, 1:2),
    ", deterministic = ",
    PCLF.is_deterministic(graph, 1:2),
)

f = rotation_system(RHO, Q)
X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
R1 = LazySets.Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
R2 = LazySets.Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])

(; quotient) = build_quotient(problem, pclf; atol = 1e-3, max_slices = 8, print_level = 0)
parts, faces = PCQ.cell_complexities(quotient)

println()
println("path-complete certificate (", Q, " nodes × 4 facets): quotient built")
@printf(
    "  %d cells, Σ pieces %d, Σ facets %d, max facets per cell %d\n",
    length(quotient.states),
    sum(parts),
    sum(faces),
    maximum(faces),
)
println("single-node certificate (1 node × 4 facets): rate ≥ 1, so no quotient exists")

display(plot_cells_by_node(quotient))
