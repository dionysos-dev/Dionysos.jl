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
using Spot
import HiGHS

const RHO = 0.94        # their published contraction rate
const GAMMA_X = 10.0    # their working set is {V ≤ Γ_X}
const GAMMA_D = 5.063   # their chosen terminal level

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
const L = gol_lazar_belta_L()

# ---------------------------------------------------------
# Their certificate, as the single-node case of our construction
# ---------------------------------------------------------

graph = PCLF.generate_DeBruijn_edges(2, 0)
pclf = PCLF.PCLF(
    graph,
    Dict{Any, PCLF.AbstractPiece}(1 => PCLF.PolyhedralPiece(L, ones(size(L, 1)))),
    RHO,
)

println(
    "graph: ",
    length(graph.verts),
    " node, complete = ",
    PCLF.is_complete(graph, 1:2),
    ", template rows = ",
    size(L, 1),
)

# ---------------------------------------------------------
# Check their certificate before using it
# ---------------------------------------------------------
# The published rate is a claim about the data; verifying it separates a misreading of the paper
# from an error in the construction that follows. The pieces are polyhedral, so the rate along an
# edge is the maximum of a linear functional over a polytope and `certify_pclf` returns the exact
# supremum rather than a sample of it.
#
# It comes out at 0.940008, marginally above the 0.94 printed in the paper: their figure is that
# supremum rounded to two decimals, not a bound the certificate meets. Nothing here depends on the
# difference — the slice recursion uses 0.94.

certified = PCLF.certify_pclf(pclf, f, HiGHS.Optimizer)

# Covering [Γ_D, Γ_X] geometrically at ratio ρ takes ⌈log(Γ_X/Γ_D)/log(1/ρ)⌉ steps.
expected_slices = ceil(Int, log(GAMMA_X / GAMMA_D) / log(1 / RHO))

println()
@printf("contraction rate: exact %.6f, paper states %.2f\n", certified.rate, RHO)
@printf("slice count:      expected %d, paper states 11\n", expected_slices)

t0 = time()
(; quotient, D) =
    build_quotient(problem, pclf; atol = 1e-3, max_slices = 40, print_level = 0)
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
# Their Problem 5.1 is synthesis: find the largest set of initial states from which *some*
# switching sequence satisfies φ. Their Problem 5.2 is verification: the set from which *every*
# switching sequence does. The two differ only in who owns the switching signal, and that is how
# they are asked here — the same quotient, specification and solver call, on the system declared
# with controlled switching for 5.1 and with autonomous switching for 5.2. Under autonomous
# switching the quotient is folded and the unchanged solver computes the universal fixed point; a
# verification returns the set and no controller.

φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)

spec = (
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
)

result = synthesize_cosafe_ltl(f, quotient, spec..., x0; N = 50)

verification = synthesize_cosafe_ltl(
    ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    quotient,
    spec...,
    x0;
    print_level = 0,
)

winning_volume =
    PCQ.get_volume(quotient, result.controllable_set; backend = CDDLib.Library())
verified_volume =
    PCQ.get_volume(quotient, verification.controllable_set; backend = CDDLib.Library())

# The soundness relation between the two quantifiers, checked rather than assumed: a state every
# switching sequence serves is in particular a state some controller wins.
@assert verification.controllable_set ⊆ result.controllable_set
@assert verification.controller === nothing

println()
@printf(
    "synthesis (5.1):    %d of %d cells winning, volume %.4f, trajectory %d steps\n",
    length(result.controllable_set),
    length(quotient.states),
    winning_volume,
    length(result.X),
)
@printf(
    "verification (5.2): %d of %d cells verified, volume %.4f\n",
    length(verification.controllable_set),
    length(quotient.states),
    verified_volume,
)

# ---------------------------------------------------------
# The gap between the quantifiers, witnessed
# ---------------------------------------------------------
# Their initial point a = (-4, -7) is controllable (5.1): some switching sequence satisfies φ
# from it. If it is not verified (5.2), a bare "no" is not evidence — the counterexample is: a
# switching word from which no controller recovers, replayed concretely.

x0_verified = any(qid -> x0 ∈ quotient.states[qid].set, verification.controllable_set)
println()
if x0_verified
    println("x0 = $x0 is verified: every switching sequence satisfies φ from it.")
    cex = nothing
else
    cex = PCQ.verification_counterexample(verification.optimizer, x0)
    tail =
        cex.entered_sink ? "runs off the abstraction's coverage" :
        "then repeats from step $(cex.lasso_start) forever"
    @printf(
        "x0 = %s is controllable but not verified. Defeating switching word: %s, %s.\n",
        string(x0),
        join(string.(cex.modes), " "),
        tail,
    )
end

# ---------------------------------------------------------
# Figures, following theirs
# ---------------------------------------------------------
# Their Figure 1 shows the working set, the terminal set and the sublevel sets; Figure 2 the
# resulting quotient; Figure 3 the set of satisfying initial states with sample trajectories.

# `show_contours` stays off: a slice is a union of polytopes, so stroking it draws the internal
# cuts of that union and the slice reads as several separate regions. The colours already separate
# consecutive slices.
fig = plot_quotient(quotient, problem, "Sublevel sets and slices"; what = :slices)
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_slices.png"))

fig = plot_quotient(quotient, problem, "Bisimulation quotient, |S| = 1")
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_quotient.png"))

# The mirrored pair the paper puts side by side: same quotient, same colour code, one word of
# difference. Synthesis (5.1, ∃) — green is what some switching sequence serves, with the
# closed-loop trajectory. Verification (5.2, ∀) — green is what every switching sequence serves,
# with the counterexample trajectory when there is one: the environment's play from the paper's
# own initial point, leaving the region no controller can bring it back from.
fig = plot_synthesis_result(quotient, result, problem, "Synthesis (5.1, ∃)")
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_synthesis.png"))

fig = plot_synthesis_result(quotient, verification, problem, "Verification (5.2, ∀)")
if cex !== nothing
    plot!(fig, ST.Trajectory([SVector{2}(x) for x in cex.X]); label = "counterexample")
end
savefig(fig, joinpath(@__DIR__, "gol_lazar_belta_verification.png"))

println("wrote gol_lazar_belta_{slices,quotient,synthesis,verification}.png")
