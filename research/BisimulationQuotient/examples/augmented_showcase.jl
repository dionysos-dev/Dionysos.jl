# The paper's worked example: a two-layer augmented quotient whose memory is provably necessary,
# and a closed-loop trajectory that hops between the layers at every step.
#
#     A₁ = ρ R(+45°),   A₂ = ρ R(−45°),   ρ = 0.80
#     graph  ℤ₂ with BOTH modes flipping the node — the node is the parity of time
#     PCLF   V₀(x) = ‖x‖_∞  on layer 0,   V₁(x) = ‖R(−45°)x‖_∞  on layer 1
#
# The certificate is period-2 and exact: R(−45°)·ρR(+45°) = ρI and R(−45°)·ρR(−45°) = ρR(−90°),
# and the box is 90°-symmetric, so every edge contracts at exactly ρ. What makes the example more
# than pretty is that the memory is *necessary at this budget*: a ±45° rotation admits no
# single-node certificate with four facets, because a convex 4-gon cannot carry the order-8
# symmetry the pair of rotations demands — the script demonstrates it by sweeping the single-node
# solver over every template orientation and certifying nothing, then certifying the two-node
# family exactly by linear programming.
#
# The two layers also *look* different — square annuli on one, diamond annuli on the other — so
# the augmented figure explains the construction by itself: one plane, two memories, and the run
# zippering between them.

include(joinpath(dirname(@__DIR__), "common.jl"))

using Spot
using Printf
import HiGHS

const ρ = 0.80
const α = π / 4

rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

f = ST.with_switching(
    HybridSystems.discreteswitchedsystem([ρ .* rotation(α), ρ .* rotation(-α)]),
    HybridSystems.ControlledSwitching(),
)

# ---------------------------------------------------------
# The period-2 certificate, and why one node cannot do it
# ---------------------------------------------------------

graph = PCLF.edgeList_to_LabDigraph([(0, 1, 1), (0, 1, 2), (1, 0, 1), (1, 0, 2)])
pclf = PCLF.PCLF(
    graph,
    Dict{Any, PCLF.AbstractPiece}(
        0 => PCLF.PolyhedralPiece(Matrix(1.0LinearAlgebra.I, 2, 2), ones(2)),
        1 => PCLF.PolyhedralPiece(rotation(-α), ones(2)),
    ),
    ρ,
)

println(
    "graph on ℤ₂: complete = ",
    PCLF.is_complete(graph, 1:2),
    ", deterministic = ",
    PCLF.is_deterministic(graph, 1:2),
)
@printf(
    "period-2 certificate: exact rate %.6f (claimed %.2f)\n",
    PCLF.certify_pclf(pclf, f, HiGHS.Optimizer).rate,
    ρ,
)

# The same budget without memory: two template rows on a single node, the orientation swept over
# a quarter turn. Nothing certifies — the necessity claim, checked rather than asserted.
let
    g0 = PCLF.generate_DeBruijn_edges(2, 0)
    sdp = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 1000,
        MOI.Silent() => true,
    )
    best = Inf
    for θ in range(0, π / 2; length = 25)[1:(end - 1)]
        p = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
            f,
            g0,
            sdp;
            Gmats = Dict(1 => rotation(θ)),
            MLF = true,
        )
        isfinite(p.JSRapprox) || continue
        best = min(best, PCLF.certify_pclf(p, f, HiGHS.Optimizer).rate)
    end
    @printf(
        "single node at the same budget: best exact rate %.4f over 24 orientations — no certificate.\n",
        best,
    )
    @assert best >= 1.0
end

# ---------------------------------------------------------
# The problem and the quotient
# ---------------------------------------------------------
# Three regions placed so that ONE tour formula has nonempty answers under BOTH quantifiers —
# every step scales the Euclidean radius by exactly ρ and turns ±45°, and each region is
# engineered against that fact (a first placement ignored the box-radius/Euclidean-radius
# factor at 45° and asked for an unrealizable spec — the solver's empty winning set was the
# reviewer):
#
#   R1, a wide sector (axis angle 0, half-width 60°, band 0.45–1.0): a visit obligation
#       survives the ∀-quantifier only where the visit is already made or forced — every
#       R1-cell satisfies `!R3 U R1` outright.
#   R2, the gateway sector just below it (same axis, half-width 60°, band 0.26–0.44): for x in
#       R1 with angle within ±15° of the axis and radius 0.46–0.78, BOTH one-step images
#       (radius 0.8 r, angles ±45°) land in R2 — |±45° + θ| ≤ 60° and the projection
#       0.8 r cos(θ ± 45°) falls in the band — so from that lens the rest of the tour is
#       forced under every word: `F(R2 & F(D))` is adversary-proof there, and the verified
#       set is exactly that lens. (A region may NOT contain the origin — the terminal level
#       τD must avoid every region — so "unavoidable by centrality" is structurally off the
#       table, and forcing through a gateway is the mechanism that remains.)
#   R3, outward (radius band 1.0–1.9 at the bottom): blocks the bottom route from x₀ (hit at
#       step 2, (0, −1.15)), so the ∃-tour must go via the top; under ∀ it is the word `1 1`
#       driving x₀ into R3 that makes x₀ controllable but NOT robust.
#
# The radius decay also makes the ∃-answer non-trivial for free: from below R2 the gateway is
# unreachable (the radius only shrinks), so the winning set has a red core on every layer.

sector(θc, w, rin, rout) = LazySets.VPolytope([
    rotation(θc) * [rin, -rin * tan(w)],
    rotation(θc) * [rin, rin * tan(w)],
    rotation(θc) * [rout, rout * tan(w)],
    rotation(θc) * [rout, -rout * tan(w)],
])

X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
R1 = sector(0.0, π / 3, 0.45, 1.00)
R2 = sector(0.0, π / 3, 0.26, 0.44)
R3 = LazySets.Hyperrectangle(; low = [-0.25, -1.9], high = [0.25, -1.0])
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])

(; quotient, D) =
    build_quotient(problem, pclf; atol = 1e-3, max_slices = 40, print_level = 1)

# ---------------------------------------------------------
# Synthesis
# ---------------------------------------------------------

φ = ltl"((!R3 U R1) & F(R2 & F(D)))"
x0 = SVector(-1.8, 0.0)

result = synthesize_cosafe_ltl(
    f,
    quotient,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    N = 40,
)

node_sequence =
    [quotient.states[m[2]].node for m in result.M if haskey(quotient.states, m[2])]
changes =
    count(i -> node_sequence[i] != node_sequence[i + 1], 1:(length(node_sequence) - 1))

parts, faces = PCQ.cell_complexities(quotient)
println()
@printf(
    "quotient: %d cells over %d nodes, Σ facets %d, max facets/cell %d\n",
    length(quotient.states),
    length(keys(quotient.slices)),
    sum(faces),
    maximum(faces)
)
@printf("trajectory: %d steps, %d node changes\n", length(result.X), changes)
println("node sequence: ", join(string.(node_sequence), " → "))

# ---------------------------------------------------------
# Verification, from the same quotient
# ---------------------------------------------------------
# One word changes the question: declaring the switching autonomous hands the modes to the
# environment, the quotient is folded, and the same solver answers what holds under EVERY
# switching sequence — for the SAME tour formula. The region placement above is what keeps this
# answer nonempty: from the gateway lens inside R1 both one-step images land in R2, so there the
# whole tour is forced under every word, and the verified set is that lens; elsewhere the
# adversary dodges R1 or detours through R3 — at x₀ the word `1 1` runs into R3 first.
# Controllable but not robust, with the certificate to show it.

verification = synthesize_cosafe_ltl(
    ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    quotient,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 1,
)
@printf(
    "verification (∀): %d of %d cells verified\n",
    length(verification.controllable_set),
    length(quotient.states),
)

x0_verified = any(qid -> x0 ∈ quotient.states[qid].set, verification.controllable_set)
cex = x0_verified ? nothing : PCQ.verification_counterexample(verification.optimizer, x0)
if cex !== nothing
    tail =
        cex.entered_sink ? "runs off the abstraction's coverage" :
        "then repeats from step $(cex.lasso_start) forever"
    @printf(
        "x0 is controllable but not verified; defeating switching word: %s, %s.\n",
        join(string.(cex.modes), " "),
        tail,
    )
end

# ---------------------------------------------------------
# Memory against geometry, at the same total budget
# ---------------------------------------------------------
# The two-node certificate spends four template rows as 2 nodes × 2 rows. The same four rows
# spent on a single node buy the octagonal gauge — and it is no ad-hoc alternative: it is the
# **common Lyapunov function the PCLF itself induces** (Angeli, Athanasopoulos, Jungers &
# Philippe, HSCC 2017), `max(V₀, V₁)`, whose sublevel set is the intersection of the square and
# the diamond and whose rows are the union of the two nodes' rows — the facet budget Σ lₛ = 4
# spent on one node instead of two. On this system the max is exactly invariant, because each
# mode swaps the two gauges: V₀(A·x) = ρV₁(x) and V₁(A·x) = ρV₀(x). Together with the failed
# 2-row sweep above, the example carries the whole trade in miniature:
#
#     1 node × 2 rows   —  certifies nothing, at any orientation
#     2 nodes × 2 rows  —  rate ρ, the augmented quotient of the figures
#     1 node × 4 rows   —  rate ρ, one layer of octagonal annuli
#
# The two successes are *different certificates*, so graph invariance does not apply — that
# theorem holds a piece family fixed — and their answers need not coincide exactly. Measured on a
# common probe grid they nearly do (about 97% agreement, small disagreements in both directions),
# while the complexity splits the other way than a "memory wins" story would like: the complete
# ℤ₂ augmentation roughly replicates the partition (≈ 2.4× the cells of the octagon), as
# `prop:complete-redundant` predicts. The example's honest moral is the paper's central one:
# at bounded per-node complexity memory buys *existence* — nothing with 2 rows certifies without
# it — not economy.

c = cos(π / 4)
G_oct = [1.0 0.0; c c; 0.0 1.0; -c c]
pclf_oct = PCLF.PCLF(
    PCLF.generate_DeBruijn_edges(2, 0),
    Dict{Any, PCLF.AbstractPiece}(1 => PCLF.PolyhedralPiece(G_oct, ones(4))),
    ρ,
)
@printf(
    "\noctagon certificate (1 node × 4 rows): exact rate %.6f\n",
    PCLF.certify_pclf(pclf_oct, f, HiGHS.Optimizer).rate,
)

# The library's own observer construction of the induced common function agrees, checked by
# sampling (its min-of-max pieces are beyond the LP certifier's polyhedral scope).
@printf(
    "induced common via build_common_lyapunov: sampled rate %.4f\n",
    PCLF.check_pclf(PCLF.build_common_lyapunov(pclf), f).rate,
)

t_oct = @elapsed oct =
    build_quotient(problem, pclf_oct; atol = 1e-3, max_slices = 40, print_level = 1)
quotient_oct = oct.quotient
D_oct = oct.D

result_oct = synthesize_cosafe_ltl(
    f,
    quotient_oct,
    Dionysos.spot_stepper(φ),
    Dict(:D => D_oct, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    N = 40,
)

# The same ∀-question asked of the memoryless construction, so the two certificates are compared
# on synthesis and verification alike.
verification_oct = synthesize_cosafe_ltl(
    ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    quotient_oct,
    Dionysos.spot_stepper(φ),
    Dict(:D => D_oct, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 1,
)
@printf(
    "common verification (∀): %d of %d cells verified\n",
    length(verification_oct.controllable_set),
    length(quotient_oct.states),
)
x0_verified_oct =
    any(qid -> x0 ∈ quotient_oct.states[qid].set, verification_oct.controllable_set)
cex_oct =
    x0_verified_oct ? nothing :
    PCQ.verification_counterexample(verification_oct.optimizer, x0)

vol(q, ids) = PCQ.get_volume(q, ids; backend = CDDLib.Library())

println()
@printf(
    "%-18s %7s %7s %9s %11s %10s %9s\n",
    "certificate",
    "cells",
    "slices",
    "Σ facets",
    "max f/cell",
    "win vol",
    "build s",
)
for (name, q, res, tb) in [
    ("2 nodes × 2 rows", quotient, result, NaN),
    ("1 node × 4 rows", quotient_oct, result_oct, t_oct),
]
    p, fc = PCQ.cell_complexities(q)
    @printf(
        "%-18s %7d %7d %9d %11d %10.3f %9.1f\n",
        name,
        length(q.states),
        length(first(values(q.slices))),
        sum(fc),
        maximum(fc),
        vol(q, res.controllable_set),
        tb,
    )
end

# The volumes above are taken over each certificate's own covered star and are therefore not
# comparable — the same trap as everywhere in this folder. The probe grid is: the same concrete
# points of the working set, asked of both winning sets.
in_win(q, res, pt) = any(qid -> pt ∈ q.states[qid].set, res.controllable_set)
probes = [
    [xx, yy] for xx in range(-1.9, 1.9; length = 39) for yy in range(-1.9, 1.9; length = 39)
]
wins_memory = count(p -> in_win(quotient, result, p), probes)
wins_geometry = count(p -> in_win(quotient_oct, result_oct, p), probes)
agree =
    count(p -> in_win(quotient, result, p) == in_win(quotient_oct, result_oct, p), probes)
@printf(
    "probe grid over X (%d points): memory wins %d, geometry wins %d, agreement %d\n",
    length(probes),
    wins_memory,
    wins_geometry,
    agree,
)

# ---------------------------------------------------------
# The figure set
# ---------------------------------------------------------
# Eight figures, one story: the problem alone, then for each certificate its ∃ and ∀ planar
# answers with the trajectory, and the 3-D views — by cell for both constructions, and the
# satisfaction per layer for the augmented one. All 2-D figures share the folder's colour
# scheme (green wins, red does not); the trajectory on a verification figure is the defeating
# counterexample, since there is no closed loop to show.

# (1) The problem by itself: the regions the specification talks about and the initial point.
# (Adding the slice family back: plot!(fig, quotient; what = :slices, fillalpha = 0.08).)
fig = plot(;
    aspect_ratio = :equal,
    legend = :outerright,
    size = (760, 520),
    title = "Regions of interest",
)
plot!(
    fig,
    problem;
    plot_region = true,
    region_alpha = 0.04,
    observation_colors = OBSERVATION_COLORS,
    observation_region_alpha = 0.3,
    observation_linewidth = 2.0,
)
scatter!(fig, [x0[1]], [x0[2]]; color = :black, markersize = 5, label = "x₀")
savefig(fig, joinpath(@__DIR__, "augmented_showcase_problem.png"))
println("wrote augmented_showcase_problem.png")

# (2)–(3) The augmented certificate's planar pair.
fig = plot_synthesis_result(quotient, result, problem, "PCLF, synthesis (∃): " * string(φ))
savefig(fig, joinpath(@__DIR__, "augmented_showcase_pclf_synthesis.png"))
println("wrote augmented_showcase_pclf_synthesis.png")

fig = plot_synthesis_result(
    quotient,
    verification,
    problem,
    "PCLF, verification (∀): " * string(φ),
)
if cex !== nothing
    plot!(fig, ST.Trajectory([SVector{2}(x) for x in cex.X]); label = "counterexample")
end
savefig(fig, joinpath(@__DIR__, "augmented_showcase_pclf_verification.png"))
println("wrote augmented_showcase_pclf_verification.png")

# (4)–(5) The induced common's planar pair.
fig = plot_synthesis_result(
    quotient_oct,
    result_oct,
    problem,
    "Common, synthesis (∃): " * string(φ),
)
savefig(fig, joinpath(@__DIR__, "augmented_showcase_common_synthesis.png"))
println("wrote augmented_showcase_common_synthesis.png")

fig = plot_synthesis_result(
    quotient_oct,
    verification_oct,
    problem,
    "Common, verification (∀): " * string(φ),
)
if cex_oct !== nothing
    plot!(fig, ST.Trajectory([SVector{2}(x) for x in cex_oct.X]); label = "counterexample")
end
savefig(fig, joinpath(@__DIR__, "augmented_showcase_common_verification.png"))
println("wrote augmented_showcase_common_verification.png")

# ---------------------------------------------------------
# Augmented 3-D views
# ---------------------------------------------------------
# Loading a Makie backend activates DionysosMakieExt; CairoMakie gives static figures for the
# paper. Everything Plots-based must come before this point — both packages export `plot`.

using CairoMakie

node_z = Dict(0 => 0.0, 1 => 1.0)

# (6) The augmented quotient with one colour per cell — the bisimulation's granularity, and the
# trajectory hopping between the layers at every step.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "memory (time parity)",
    zticks = ([0.0, 1.0], ["0", "1"]),
    title = "Augmented quotient, one colour per cell",
)
DI.plot_augmented_bisimulation!(
    ax,
    quotient;
    node_z = node_z,
    color_by = :state,
    alpha = 0.35,
    show_contours = false,
)
DI.plot_augmented_trajectory!(ax, quotient, result.X, result.M; node_z = node_z)
CairoMakie.save(joinpath(@__DIR__, "augmented_showcase_pclf_3d_cells.png"), mk)
println("wrote augmented_showcase_pclf_3d_cells.png")

# (7) The same by-cell view of the memoryless construction — one flat layer, the run never
# leaves the plane: what the augmentation adds is literally the vertical axis.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "no memory",
    zticks = ([0.0, 1.0], ["", ""]),
    title = "Induced common (one node), one colour per cell",
)
node_z_oct = Dict(1 => 0.0)
DI.plot_augmented_bisimulation!(
    ax,
    quotient_oct;
    node_z = node_z_oct,
    color_by = :state,
    alpha = 0.35,
    show_contours = false,
)
DI.plot_augmented_trajectory!(
    ax,
    quotient_oct,
    result_oct.X,
    result_oct.M;
    node_z = node_z_oct,
)
CairoMakie.save(joinpath(@__DIR__, "augmented_showcase_common_3d_cells.png"), mk)
println("wrote augmented_showcase_common_3d_cells.png")

# (8) The satisfaction itself, augmented: winning cells green, losing red, per layer. This is
# the figure the planar view cannot draw — the winning set lives on *augmented* states, so the
# same region of the plane can be winning on one memory and losing on the other, and here is
# where.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "memory (time parity)",
    zticks = ([0.0, 1.0], ["0", "1"]),
    title = "Satisfaction per layer (∃): " * string(φ),
)
DI.plot_augmented_bisimulation!(
    ax,
    quotient;
    state_ids = result.uncontrollable_set,
    node_z = node_z,
    color_by = LOSING_COLOR,
    alpha = 0.30,
    show_contours = false,
)
DI.plot_augmented_bisimulation!(
    ax,
    quotient;
    state_ids = result.controllable_set,
    node_z = node_z,
    color_by = WINNING_COLOR,
    alpha = 0.45,
    show_contours = false,
)
DI.plot_augmented_trajectory!(ax, quotient, result.X, result.M; node_z = node_z)
CairoMakie.save(joinpath(@__DIR__, "augmented_showcase_pclf_3d_satisfaction.png"), mk)
println("wrote augmented_showcase_pclf_3d_satisfaction.png")

# (9) The satisfaction of the VERIFICATION question, per layer: green survives every switching
# word, red does not. Nonempty green is what makes the panel worth drawing, and the
# counterexample from x0 threads the layers like the closed loop does, alternating parity.
function qid_at(x, node)
    # findfirst returns a position, and the state ids are not dense - return the id itself.
    ks = sort(collect(keys(quotient.states)))
    i = findfirst(
        qid -> quotient.states[qid].node == node && x in quotient.states[qid].set,
        ks,
    )
    return i === nothing ? nothing : ks[i]
end

mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x1",
    ylabel = "x2",
    zlabel = "memory (time parity)",
    zticks = ([0.0, 1.0], ["0", "1"]),
    title = "Satisfaction per layer (forall): " * string(φ),
)
DI.plot_augmented_bisimulation!(
    ax,
    quotient;
    state_ids = verification.uncontrollable_set,
    node_z = node_z,
    color_by = LOSING_COLOR,
    alpha = 0.30,
    show_contours = false,
)
DI.plot_augmented_bisimulation!(
    ax,
    quotient;
    state_ids = verification.controllable_set,
    node_z = node_z,
    color_by = WINNING_COLOR,
    alpha = 0.45,
    show_contours = false,
)
if cex !== nothing
    X_cex = SVector{2, Float64}[]
    M_cex = Tuple{Int, Int}[]
    for (i, x) in enumerate(cex.X)
        qid = qid_at(SVector{2}(x), mod(i - 1, 2))
        qid === nothing && break
        push!(X_cex, SVector{2}(x))
        push!(M_cex, (0, qid))
    end
    isempty(X_cex) ||
        DI.plot_augmented_trajectory!(ax, quotient, X_cex, M_cex; node_z = node_z)
end
CairoMakie.save(
    joinpath(@__DIR__, "augmented_showcase_pclf_3d_satisfaction_forall.png"),
    mk,
)
println("wrote augmented_showcase_pclf_3d_satisfaction_forall.png")
