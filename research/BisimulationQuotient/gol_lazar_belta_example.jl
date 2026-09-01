# Reproduction of Example 3.1 of Gol, Ding, Lazar & Belta, *Finite Bisimulations for Switched
# Linear Systems* (arXiv:1208.5471).
#
# Their construction is the |S| = 1 case of the one developed here: with a single node the
# path-complete graph is one node carrying a self-loop per mode, the node-dependent slice family
# collapses to a single family, and the refinement reduces to their Algorithm 2. Running their
# published certificate through this implementation should return their published numbers, which
# is what makes the claim that the framework generalises theirs checkable rather than rhetorical.
#
# Their data, from Examples 3.1 and 4.1:
#     A₁ = [-0.65 0.32; -0.42 -0.92],  A₂ = [0.65 0.32; -0.42 -0.92]
#     V(x) = ‖Lx‖_∞ with the four-row L below, contraction rate ρ = 0.94
#     Γ_X = 10, Γ_D = 5.063, and consequently N = 11 slices.
#
# The certificate is theirs and is written down rather than recomputed, so what is being tested is
# our construction and not a solver.
#
# This script deliberately does **not** compare their certificate against a path-complete one. A
# comparison needs both quotients to partition the same region at comparable resolution, and on
# this example the conic templates we would otherwise use produce a terminal set that absorbs most
# of the domain — a single-node conic-order-1 certificate yields a quotient of one cell. Cell
# counts across such certificates measure the terminal set, not the abstraction. The regime where
# memory does pay is exhibited by `memory_vs_geometry.jl` instead.

include(joinpath(@__DIR__, "common.jl"))

using Printf
using Random

const RHO = 0.94        # their published contraction rate
const GAMMA_X = 10.0    # their working set is {V ≤ Γ_X}
const GAMMA_D = 5.063   # their chosen terminal level

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
const L = gol_lazar_belta_L()

# ---------------------------------------------------------
# Check their certificate before using it
# ---------------------------------------------------------
# The published rate is a claim about the data; verifying it separates a misreading of the paper
# from an error in the construction that follows.

function measured_rate(; nsample = 50_000, seed = 1)
    modes = ST.mode_matrices(f)
    V(x) = maximum(abs.(L * x))
    rng = Random.MersenneTwister(seed)
    worst = 0.0
    for _ in 1:nsample
        x = randn(rng, 2)
        v = V(x)
        v < 1e-12 && continue
        worst = max(worst, maximum(V(A * x) for A in modes) / v)
    end
    return worst
end

# Covering [Γ_D, Γ_X] geometrically at ratio ρ takes ⌈log(Γ_X/Γ_D)/log(1/ρ)⌉ steps.
expected_slices = ceil(Int, log(GAMMA_X / GAMMA_D) / log(1 / RHO))

@printf("contraction rate: measured %.4f, paper states %.2f\n", measured_rate(), RHO)
@printf("slice count:      expected %d, paper states 11\n", expected_slices)

# ---------------------------------------------------------
# Their certificate, as the single-node case of our construction
# ---------------------------------------------------------

graph = PCLF.generate_DeBruijn_edges(2, 0)
pclf = PCLF.PCLF(
    graph,
    Dict{Any, PCLF.AbstractPiece}(1 => PCLF.PolyhedralPiece(L, ones(size(L, 1)))),
    RHO,
)

println()
println(
    "graph: ",
    length(graph.verts),
    " node, complete = ",
    PCLF.is_complete(graph, 1:2),
    ", template rows = ",
    size(L, 1),
)

t0 = time()
(; quotient) = build_quotient(problem, pclf; atol = 1e-3, max_slices = 40, print_level = 0)
build_s = time() - t0

parts, faces = PCQ.cell_complexities(quotient)
n_slices = length(first(values(quotient.slices)))

println()
@printf("slices built: %d\n", n_slices)
@printf(
    "quotient: %d cells over %d node, Σ pieces %d, Σ facets %d, max facets/cell %d, %.1f s\n",
    length(quotient.states),
    length(keys(quotient.slices)),
    sum(parts),
    sum(faces),
    maximum(faces),
    build_s,
)

# ---------------------------------------------------------
# Figures, following theirs
# ---------------------------------------------------------
# Their Figure 1 shows the working set, the terminal set and the sublevel sets; their Figure 2 the
# resulting quotient.

fig = plot(; aspect_ratio = :equal, legend = false, title = "Sublevel sets and slices")
plot!(fig, quotient; what = :slices, show_contours = true)
plot!(fig, problem; opacity = 0.25)
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_slices.png"))

fig =
    plot(; aspect_ratio = :equal, legend = false, title = "Bisimulation quotient, |S| = 1")
plot!(fig, quotient; what = :states, show_contours = false)
plot!(fig, problem; opacity = 0.25)
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_quotient.png"))

println("wrote gol_lazar_belta_slices.png and gol_lazar_belta_quotient.png")
