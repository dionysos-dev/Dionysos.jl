# The second worked example: an *incomplete* path-complete graph that beats its own induced
# common function on every complexity count, at identical pieces and identical certified rate.
#
#     A₁ = ρ R(+45°),   A₂ = ρ R(−45°),   ρ = 0.80        — the system of `augmented_showcase.jl`
#     graph  dual De Bruijn of order 1: node (i) plays only mode i, then moves to any node
#     pieces the SAME octagonal gauge on both nodes
#
# The graph is incomplete — node (i) has no edge for the other mode — but path-complete: a word
# is followed by standing, at every step, on the node that announces it. The node therefore means
# "the next mode", and the augmented trajectory jumps between layers exactly at the *decision
# points* where the controller changes mode — where the first example's time-parity graph jumped
# at every step.
#
# Because the two pieces are identical, the induced common Lyapunov function is that piece
# itself, and the single-node run *is* the induced-common construction: the comparison below
# isolates the pure effect of the graph at fixed gauge, fixed rate and fixed slice family. A
# complete graph would replicate the partition |S| times (`prop:complete-redundant`); the
# incomplete one prunes each layer to its own mode's transitions instead, and the quotient comes
# out *smaller* than the common's — the reduction first observed on the De Bruijn duals and so
# far unexplained, here put to work. With identical pieces the synthesis answers must agree
# (the invariance theorems' regime), which the probe grid checks.

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

c = cos(π / 4)
G_oct = [1.0 0.0; c c; 0.0 1.0; -c c]
octagon() = PCLF.PolyhedralPiece(G_oct, ones(4))

# ---------------------------------------------------------
# The two constructions: one node, and the incomplete augmentation of the same gauge
# ---------------------------------------------------------

pclf_common = PCLF.PCLF(
    PCLF.generate_DeBruijn_edges(2, 0),
    Dict{Any, PCLF.AbstractPiece}(1 => octagon()),
    ρ,
)

dual_graph = PCLF.generate_DeBruijn_edges(2, 1; dual = true)
pclf_dual = PCLF.PCLF(
    dual_graph,
    Dict{Any, PCLF.AbstractPiece}((1,) => octagon(), (2,) => octagon()),
    ρ,
)

println(
    "dual De Bruijn: complete = ",
    PCLF.is_complete(dual_graph, 1:2),
    ", deterministic = ",
    PCLF.is_deterministic(dual_graph, 1:2),
    ", edges = ",
    length(dual_graph.edges),
)
@printf(
    "induced common (1 node): exact rate %.6f;  incomplete augmentation (2 nodes): exact rate %.6f\n",
    PCLF.certify_pclf(pclf_common, f, HiGHS.Optimizer).rate,
    PCLF.certify_pclf(pclf_dual, f, HiGHS.Optimizer).rate,
)

# ---------------------------------------------------------
# Problem: the regions of the first example, with R2 pushed outward
# ---------------------------------------------------------

X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
R1 = LazySets.Hyperrectangle(; low = [0.42, -0.16], high = [0.75, 0.16])
R2 = LazySets.Hyperrectangle(; low = [-0.16, 0.30], high = [0.16, 0.50])
R3 = LazySets.Hyperrectangle(; low = [-0.2, -1.4], high = [0.2, -0.5])
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])

build(pclf) = begin
    t = @elapsed r =
        build_quotient(problem, pclf; atol = 1e-3, max_slices = 40, print_level = 0)
    (; r.quotient, r.D, t)
end

# A coarse throwaway build first, so neither timed build pays Julia's compilation — the first
# quotient of a session otherwise reports compile time as construction time.
build_quotient(problem, pclf_common; atol = 1e-2, max_slices = 3, print_level = 0)

common = build(pclf_common)
dual = build(pclf_dual)

# ---------------------------------------------------------
# Synthesis on both, same specification
# ---------------------------------------------------------

φ = ltl"((!R3 U R1) & F(R2 & F(D)))"
x0 = SVector(-1.4, 0.0)

solve(q, D) = synthesize_cosafe_ltl(
    f,
    q,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    N = 40,
    print_level = 0,
)

result_common = solve(common.quotient, common.D)
result_dual = solve(dual.quotient, dual.D)

node_sequence = [
    dual.quotient.states[m[2]].node for
    m in result_dual.M if haskey(dual.quotient.states, m[2])
]
changes =
    count(i -> node_sequence[i] != node_sequence[i + 1], 1:(length(node_sequence) - 1))
@printf(
    "\ntrajectory on the augmentation: %d steps, %d node changes (a jump = a mode decision)\n",
    length(result_dual.X),
    changes,
)
println("announced-mode sequence: ", join(string.(first.(node_sequence)), " → "))

# ---------------------------------------------------------
# The comparison, and the agreement check
# ---------------------------------------------------------

println()
@printf(
    "%-24s %7s %7s %9s %11s %9s\n",
    "construction",
    "cells",
    "slices",
    "Σ facets",
    "max f/cell",
    "build s",
)
for (name, b) in [("induced common (1 node)", common), ("incomplete PCLF (2)", dual)]
    p, fc = PCQ.cell_complexities(b.quotient)
    @printf(
        "%-24s %7d %7d %9d %11d %9.1f\n",
        name,
        length(b.quotient.states),
        length(first(values(b.quotient.slices))),
        sum(fc),
        maximum(fc),
        b.t,
    )
end

# Identical pieces put the two runs in the invariance theorems' regime: the certified answers
# must coincide, not merely resemble each other.
in_win(q, res, pt) = any(qid -> pt ∈ q.states[qid].set, res.controllable_set)
probes = [
    [xx, yy] for xx in range(-1.9, 1.9; length = 39) for yy in range(-1.9, 1.9; length = 39)
]
agree = count(
    p ->
        in_win(common.quotient, result_common, p) == in_win(dual.quotient, result_dual, p),
    probes,
)
@printf(
    "probe grid over X (%d points): common wins %d, augmentation wins %d, agreement %d\n",
    length(probes),
    count(p -> in_win(common.quotient, result_common, p), probes),
    count(p -> in_win(dual.quotient, result_dual, p), probes),
    agree,
)

# ---------------------------------------------------------
# Verification (∀) on both constructions
# ---------------------------------------------------------
# The switching handed to the environment: what holds under EVERY switching word. On this
# system the adversary can spin the state away from R1 forever, so both answers should be
# near-empty, with a lasso counterexample from x₀ naming the defeating word.

verify(q, D) = synthesize_cosafe_ltl(
    ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    q,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)

verification_common = verify(common.quotient, common.D)
verification_dual = verify(dual.quotient, dual.D)

counterexample(q, verification) =
    any(qid -> x0 ∈ q.states[qid].set, verification.controllable_set) ? nothing :
    PCQ.verification_counterexample(verification.optimizer, x0)

cex_common = counterexample(common.quotient, verification_common)
cex_dual = counterexample(dual.quotient, verification_dual)

@printf(
    "verification (∀): common %d of %d cells, augmentation %d of %d cells\n",
    length(verification_common.controllable_set),
    length(common.quotient.states),
    length(verification_dual.controllable_set),
    length(dual.quotient.states),
)
if cex_dual !== nothing
    @printf(
        "defeating switching word from x₀ (augmentation): %s\n",
        join(string.(cex_dual.modes), " "),
    )
end

# ---------------------------------------------------------
# The figure set
# ---------------------------------------------------------
# The same eight figures as `augmented_showcase.jl`, same names after the prefix, so the two
# examples drop into the same layout: the problem alone, each construction's ∃/∀ planar pair
# with its trajectory (closed loop for ∃, defeating counterexample for ∀), the by-cell 3-D
# views, and the satisfaction per announcement layer.

# (1) The problem by itself.
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
savefig(fig, joinpath(@__DIR__, "incomplete_showcase_problem.png"))
println("wrote incomplete_showcase_problem.png")

# (2)–(3) The augmentation's planar pair.
fig = plot_synthesis_result(
    dual.quotient,
    result_dual,
    problem,
    "PCLF, synthesis (∃): " * string(φ),
)
savefig(fig, joinpath(@__DIR__, "incomplete_showcase_pclf_synthesis.png"))
println("wrote incomplete_showcase_pclf_synthesis.png")

fig = plot_synthesis_result(
    dual.quotient,
    verification_dual,
    problem,
    "PCLF, verification (∀): " * string(φ),
)
if cex_dual !== nothing
    plot!(fig, ST.Trajectory([SVector{2}(x) for x in cex_dual.X]); label = "counterexample")
end
savefig(fig, joinpath(@__DIR__, "incomplete_showcase_pclf_verification.png"))
println("wrote incomplete_showcase_pclf_verification.png")

# (4)–(5) The induced common's planar pair.
fig = plot_synthesis_result(
    common.quotient,
    result_common,
    problem,
    "Common, synthesis (∃): " * string(φ),
)
savefig(fig, joinpath(@__DIR__, "incomplete_showcase_common_synthesis.png"))
println("wrote incomplete_showcase_common_synthesis.png")

fig = plot_synthesis_result(
    common.quotient,
    verification_common,
    problem,
    "Common, verification (∀): " * string(φ),
)
if cex_common !== nothing
    plot!(
        fig,
        ST.Trajectory([SVector{2}(x) for x in cex_common.X]);
        label = "counterexample",
    )
end
savefig(fig, joinpath(@__DIR__, "incomplete_showcase_common_verification.png"))
println("wrote incomplete_showcase_common_verification.png")

# ---------------------------------------------------------
# Augmented 3-D views
# ---------------------------------------------------------
# Everything Plots-based must come before this point — both packages export `plot`.

using CairoMakie

node_z = Dict((1,) => 0.0, (2,) => 1.0)

# (6) The augmentation with one colour per cell — the granularity of the two announcement
# layers, and the trajectory jumping exactly at the mode decisions.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "announced next mode",
    zticks = ([0.0, 1.0], ["mode 1", "mode 2"]),
    title = "Incomplete augmentation, one colour per cell",
)
DI.plot_augmented_bisimulation!(
    ax,
    dual.quotient;
    node_z = node_z,
    color_by = :state,
    alpha = 0.35,
    show_contours = false,
)
DI.plot_augmented_trajectory!(
    ax,
    dual.quotient,
    result_dual.X,
    result_dual.M;
    node_z = node_z,
)
CairoMakie.save(joinpath(@__DIR__, "incomplete_showcase_pclf_3d_cells.png"), mk)
println("wrote incomplete_showcase_pclf_3d_cells.png")

# (7) The same by-cell view of the memoryless construction — one flat layer.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "no memory",
    zticks = ([0.0, 1.0], ["", ""]),
    title = "Induced common (one node), one colour per cell",
)
node_z_common = Dict(1 => 0.0)
DI.plot_augmented_bisimulation!(
    ax,
    common.quotient;
    node_z = node_z_common,
    color_by = :state,
    alpha = 0.35,
    show_contours = false,
)
DI.plot_augmented_trajectory!(
    ax,
    common.quotient,
    result_common.X,
    result_common.M;
    node_z = node_z_common,
)
CairoMakie.save(joinpath(@__DIR__, "incomplete_showcase_common_3d_cells.png"), mk)
println("wrote incomplete_showcase_common_3d_cells.png")

# (8) Satisfaction per announcement layer: winning green, losing red — where committing to a
# next mode decides the specification.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "announced next mode",
    zticks = ([0.0, 1.0], ["mode 1", "mode 2"]),
    title = "Satisfaction per layer (∃): " * string(φ),
)
DI.plot_augmented_bisimulation!(
    ax,
    dual.quotient;
    state_ids = result_dual.uncontrollable_set,
    node_z = node_z,
    color_by = LOSING_COLOR,
    alpha = 0.30,
    show_contours = false,
)
DI.plot_augmented_bisimulation!(
    ax,
    dual.quotient;
    state_ids = result_dual.controllable_set,
    node_z = node_z,
    color_by = WINNING_COLOR,
    alpha = 0.45,
    show_contours = false,
)
DI.plot_augmented_trajectory!(
    ax,
    dual.quotient,
    result_dual.X,
    result_dual.M;
    node_z = node_z,
)
CairoMakie.save(joinpath(@__DIR__, "incomplete_showcase_pclf_3d_satisfaction.png"), mk)
println("wrote incomplete_showcase_pclf_3d_satisfaction.png")
