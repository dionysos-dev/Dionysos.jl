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
using Spot

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
(; quotient, D) = build_quotient(problem, pclf; atol = 1e-3, max_slices = 40, print_level = 0)
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
# Their control problem
# ---------------------------------------------------------
# Their Example 5.1 states the specification in words — "a trajectory never visits R₂ and
# eventually visits R₁; moreover, if it visits R₃ then it must not visit R₁ at the next step" —
# and formalises it as
#
#     φ := (¬R₂ U Π_D) ∧ F R₁ ∧ ((R₃ ⇒ X ¬R₁) U Π_D)
#
# Their Problem 5.1 is synthesis: find the largest set of initial states from which some switching
# sequence satisfies φ. That is the existential semantics, and it is what `synthesize_cosafe_ltl`
# computes here.
#
# Their Problem 5.2 — verification under arbitrary switching, where *every* switching sequence must
# satisfy φ — is a different quantifier, over switching sequences rather than over graph nodes, and
# is not exposed by the optimizer used here. It is left out rather than approximated.

φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)

result = synthesize_cosafe_ltl(
    f,
    quotient,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    N = 50,
)

winning_volume =
    PCQ.get_volume(quotient, result.controllable_set; backend = CDDLib.Library())

println()
@printf(
    "synthesis: %d of %d cells winning, volume %.4f, trajectory %d steps\n",
    length(result.controllable_set), length(quotient.states),
    winning_volume, length(result.X),
)

# ---------------------------------------------------------
# Figures, following theirs
# ---------------------------------------------------------
# Their Figure 1 shows the working set, the terminal set and the sublevel sets; Figure 2 the
# resulting quotient; Figure 3 the set of satisfying initial states with sample trajectories.

fig = plot(; aspect_ratio = :equal, legend = false, title = "Sublevel sets and slices")
plot!(fig, quotient; what = :slices, show_contours = true)
plot!(fig, problem; opacity = 0.25)
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_slices.png"))

fig =
    plot(; aspect_ratio = :equal, legend = false, title = "Bisimulation quotient, |S| = 1")
plot!(fig, quotient; what = :states, show_contours = false)
plot!(fig, problem; opacity = 0.25)
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_quotient.png"))

fig = plot(; aspect_ratio = :equal, legend = false, title = "Satisfying initial states")
plot!(
    fig,
    quotient;
    what = :states,
    state_ids = result.controllable_set,
    show_contours = false,
    user_color = :mediumpurple,
    fillalpha = 0.95,
)
plot!(fig, problem; region_alpha = 0.0, observation_region_alpha = 0.35, plot_region = false)
plot!(fig, ST.Trajectory(result.X); label = "trajectory")
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_satisfying.png"))

println("wrote gol_lazar_belta_{slices,quotient,satisfying}.png")
