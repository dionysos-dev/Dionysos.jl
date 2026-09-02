# Example 3.1 — the paper's reference run.
#
# The system, working set and observation regions are those of Gol, Ding, Lazar & Belta
# (arXiv:1208.5471, Example 3.1); the certificate is ours. Where
# `gol_lazar_belta_example.jl` reproduces *their* certificate at |S| = 1 to check that the
# framework contains theirs, this script runs the full pipeline on the same problem with a
# two-node path-complete certificate: quotient, co-safe LTL synthesis, controller, closed-loop
# trajectory and winning-region volume.
#
# Reference numbers: 10611 quotient states, controllable-set volume 186.72811640312693.
#
# The two are not comparable as a cost benchmark. They use different certificates, and a
# certificate determines its own slice family and terminal set, so their cell counts measure
# different partitions of different regions — see the note at the head of
# `gol_lazar_belta_example.jl`.

include(joinpath(dirname(@__DIR__), "common.jl"))
using Spot

gr()

# ---------------------------------------------------------
# System and problem
# ---------------------------------------------------------

# The system, working set and observation regions are Example 3.1 of Gol, Ding, Lazar & Belta
# (arXiv:1208.5471). They are defined once in `common.jl` so that this script and
# `gol_lazar_belta_example.jl`, which compares our certificate against theirs, abstract exactly
# the same problem — otherwise the two quotients would partition different sets.

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()

# ---------------------------------------------------------
# Polyhedral PCLF
# ---------------------------------------------------------
# De Bruijn graph of order 1 over 2 modes, so two nodes, each carrying a piece built from the
# order-2 conic partition of the plane.

graph = PCLF.generate_DeBruijn_edges(2, 1)
println(graph)

partition = PCLF.conic_partitions_2d(2)
pclf = PCLF.compute_polyhedral_pieces_pclf(
    f,
    graph,
    JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000),
    Dict((1,) => partition, (2,) => partition);
    MLF = true,
)
println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

# ---------------------------------------------------------
# Quotient
# ---------------------------------------------------------
# `atol` is the inset applied when one polytope is cut from another. At 1e-3 it silently drops
# ~3.7% of the domain; at 1e-6 flat pieces survive and break plotting. 1e-4 loses ~0.4%.

(; optimizer, quotient, D) =
    build_quotient(problem, pclf; atol = 1e-4, level_tol = 1e-2, max_slices = 50)

export_optimizer_jld2(optimizer, joinpath(@__DIR__, "gol_lazar_belta_pclf.jld2"))

using Printf

# ---------------------------------------------------------
# Co-safe LTL synthesis (∃) and verification (∀)
# ---------------------------------------------------------

φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"
x0 = SVector(-4.0, -7.0)   # initial point a in the paper

result = synthesize_cosafe_ltl(
    f,
    quotient,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0,
)

println(
    "Volume of controllable set = ",
    PCQ.get_volume(quotient, result.controllable_set; backend = CDDLib.Library()),
)

verification = synthesize_cosafe_ltl(
    ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    quotient,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)

# ---------------------------------------------------------
# The memoryless arm: their own published certificate
# ---------------------------------------------------------
# The natural |S| = 1 comparison is not a solver-found single-node certificate — a conic
# single-node one collapses the quotient to almost nothing (see the README) — but the paper's
# own 4-row certificate, the one `gol_lazar_belta_example.jl` validates against their numbers.

pclf_common = PCLF.PCLF(
    PCLF.generate_DeBruijn_edges(2, 0),
    Dict{Any, PCLF.AbstractPiece}(
        1 => PCLF.PolyhedralPiece(gol_lazar_belta_L(), ones(size(gol_lazar_belta_L(), 1))),
    ),
    0.94,
)

common = build_quotient(problem, pclf_common; atol = 1e-3, max_slices = 40, print_level = 0)
quotient_common = common.quotient
D_common = common.D

result_common = synthesize_cosafe_ltl(
    f,
    quotient_common,
    Dionysos.spot_stepper(φ),
    Dict(:D => D_common, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)
verification_common = synthesize_cosafe_ltl(
    ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    quotient_common,
    Dionysos.spot_stepper(φ),
    Dict(:D => D_common, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)

@printf(
    "PCLF:   %d cells, ∃ wins %d, ∀ verifies %d\n",
    length(quotient.states),
    length(result.controllable_set),
    length(verification.controllable_set),
)
@printf(
    "common: %d cells, ∃ wins %d, ∀ verifies %d\n",
    length(quotient_common.states),
    length(result_common.controllable_set),
    length(verification_common.controllable_set),
)

# The trajectory a verification figure carries: from a verified x₀ every switching word
# satisfies φ, so an arbitrary one — the alternating word — is a witness run; from an
# unverified x₀ the defeating word of the counterexample is shown instead.
is_verified(q, verif) = any(qid -> x0 ∈ q.states[qid].set, verif.controllable_set)

function verification_trajectory(q, verif, D_set)
    if is_verified(q, verif)
        A = UT.mode_matrices(f)
        xs = [x0]
        for k in 1:30
            last(xs) ∈ D_set && break
            push!(xs, SVector{2}(A[mod1(k, length(A))] * last(xs)))
        end
        return (xs, "witness run (alternating word)")
    end
    cex = PCQ.verification_counterexample(verif.optimizer, x0)
    return ([SVector{2}(x) for x in cex.X], "counterexample")
end

# ---------------------------------------------------------
# The figure set
# ---------------------------------------------------------
# The poster figures, prefix `gol_lazar_belta_pclf_`. Every figure is written as PNG and as
# PDF — the poster wants vector graphics.

save_fig(fig, name) = begin
    savefig(fig, joinpath(@__DIR__, name * ".png"))
    savefig(fig, joinpath(@__DIR__, name * ".pdf"))
    return println("wrote ", name, ".{png,pdf}")
end

_centroid(R) = sum(LazySets.vertices_list(R)) / length(LazySets.vertices_list(R))

# The problem figure needs the induced CLF's largest slice as its backdrop, so it is drawn
# after the CLF arm below.

# (1b) The same scene over the common certificate's slice family, coloured by slice — the
# onion the regions cut, which is where the quotient's cells come from.
fig = plot(; aspect_ratio = :equal, legend = false, size = (700, 620), title = "")
plot!(fig, quotient_common; what = :slices, show_contours = false, fillalpha = 0.35)
plot!(
    fig,
    problem;
    plot_region = false,
    observation_colors = OBSERVATION_COLORS,
    observation_region_alpha = 0.3,
    observation_linewidth = 2.5,
)
for P in LazySets.array(D_common)
    plot!(fig, P; fillalpha = 0.0, linecolor = :black, linestyle = :dash, linewidth = 1.5)
end
scatter!(fig, [x0[1]], [x0[2]]; color = :black, markersize = 5)
for (name, R, c) in [
    ("R1", R1, OBSERVATION_COLORS[1]),
    ("R2", R2, OBSERVATION_COLORS[2]),
    ("R3", R3, OBSERVATION_COLORS[3]),
]
    ctr = _centroid(R)
    annotate!(fig, ctr[1], ctr[2], Plots.text(name, 16, c, :center))
end
Dctr_c = _centroid(LazySets.array(D_common)[1])
annotate!(fig, Dctr_c[1], Dctr_c[2], Plots.text("D", 16, :black, :center))
annotate!(fig, x0[1] - 0.6, x0[2] - 0.7, Plots.text("x₀", 15, :black, :center))
save_fig(fig, "gol_lazar_belta_pclf_common_slices")

# (2)–(3) The PCLF's planar pair.
fig = plot_synthesis_result(quotient, result, problem, "PCLF, synthesis (∃): " * string(φ))
save_fig(fig, "gol_lazar_belta_pclf_pclf_synthesis")

fig = plot_synthesis_result(
    quotient,
    verification,
    problem,
    "PCLF, verification (∀): " * string(φ),
)
traj, lbl = verification_trajectory(quotient, verification, D)
plot!(fig, ST.Trajectory(traj); label = lbl)
save_fig(fig, "gol_lazar_belta_pclf_pclf_verification")

# (4)–(5) Their certificate's planar pair.
fig = plot_synthesis_result(
    quotient_common,
    result_common,
    problem,
    "Common, synthesis (∃): " * string(φ),
)
save_fig(fig, "gol_lazar_belta_pclf_common_synthesis")

fig = plot_synthesis_result(
    quotient_common,
    verification_common,
    problem,
    "Common, verification (∀): " * string(φ),
)
traj, lbl = verification_trajectory(quotient_common, verification_common, D_common)
plot!(fig, ST.Trajectory(traj); label = lbl)
save_fig(fig, "gol_lazar_belta_pclf_common_verification")

# ---------------------------------------------------------
# The induced common CLF: where the complexity concentrates
# ---------------------------------------------------------
# The third arm is the common Lyapunov function the PCLF itself induces
# (`build_common_lyapunov`, the observer construction of Angeli–Athanasopoulos–Jungers–
# Philippe): the same certified information with the memory folded away into a single, much
# more complicated piece. The quotient it builds carries the point the other two arms cannot
# show — the facet budget is not removed by memory, it is *redistributed*: the CLF quotient
# concentrates it into few fat cells where the PCLF spreads it thinly.

clf = PCLF.build_common_lyapunov(pclf)
common_clf = build_quotient(problem, clf; atol = 1e-4, level_tol = 1e-2, max_slices = 50)
export_optimizer_jld2(common_clf.optimizer, joinpath(@__DIR__, "gol_lazar_belta_clf.jld2"))
quotient_clf = common_clf.quotient
D_clf = common_clf.D

result_clf = synthesize_cosafe_ltl(
    f,
    quotient_clf,
    Dionysos.spot_stepper(φ),
    Dict(:D => D_clf, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)
verification_clf = synthesize_cosafe_ltl(
    ST.with_switching(f, HybridSystems.AutonomousSwitching()),
    quotient_clf,
    Dionysos.spot_stepper(φ),
    Dict(:D => D_clf, :R1 => R1, :R2 => R2, :R3 => R3),
    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
    x0;
    print_level = 0,
)

println()
@printf(
    "%-22s %7s %10s %12s %8s %8s\n",
    "certificate",
    "cells",
    "Σ facets",
    "max f/cell",
    "∃ wins",
    "∀ ok",
)
for (name, q, res, verif) in [
    ("PCLF (2 nodes)", quotient, result, verification),
    ("their CLF (4 rows)", quotient_common, result_common, verification_common),
    ("induced CLF", quotient_clf, result_clf, verification_clf),
]
    p, fc = PCQ.cell_complexities(q)
    @printf(
        "%-22s %7d %10d %12d %8d %8d\n",
        name,
        length(q.states),
        sum(fc),
        maximum(fc),
        length(res.controllable_set),
        length(verif.controllable_set),
    )
end

# The E8 histogram: per-cell complexity of the PCLF against the common function it induces —
# the same facet budget spent thin against spent concentrated. Log counts, because most cells
# are simple and a long tail is not, and the tail is the point.
parts_pclf, facets_pclf = PCQ.cell_complexities(quotient)
parts_clf, facets_clf = PCQ.cell_complexities(quotient_clf)

function complexity_panel(pclf_values, clf_values; xlabel, title)
    bins = range(0, 1.02 * maximum(vcat(pclf_values, clf_values)); length = 55)
    fig = histogram(
        pclf_values;
        bins = bins,
        yscale = :log10,
        alpha = 0.55,
        label = "PCLF (2 nodes)",
        color = :steelblue,
        linecolor = :steelblue,
        xlabel = xlabel,
        ylabel = "number of cells",
        title = title,
        ylims = (0.7, 2e4),
        legend = :topright,
        framestyle = :box,
        grid = false,
    )
    histogram!(
        fig,
        clf_values;
        bins = bins,
        yscale = :log10,
        alpha = 0.55,
        label = "induced CLF (1 node)",
        color = :firebrick,
        linecolor = :firebrick,
    )
    # The two maxima, which is where the separation actually lives.
    vline!(
        fig,
        [maximum(pclf_values)];
        color = :steelblue,
        ls = :dash,
        lw = 1.5,
        label = "",
    )
    vline!(fig, [maximum(clf_values)]; color = :firebrick, ls = :dash, lw = 1.5, label = "")
    return fig
end

fig = plot(
    complexity_panel(
        parts_pclf,
        parts_clf;
        xlabel = "polytopes per cell",
        title = "Pieces per semilinear cell",
    ),
    complexity_panel(
        facets_pclf,
        facets_clf;
        xlabel = "facets per cell",
        title = "Facets per semilinear cell",
    );
    layout = (1, 2),
    size = (900, 340),
    left_margin = 5Plots.mm,
    bottom_margin = 5Plots.mm,
    top_margin = 3Plots.mm,
)
save_fig(fig, "gol_lazar_belta_pclf_complexity_histogram")

# (10)–(11) The induced CLF's planar pair.
fig = plot_synthesis_result(
    quotient_clf,
    result_clf,
    problem,
    "Induced CLF, synthesis (∃): " * string(φ),
)
save_fig(fig, "gol_lazar_belta_pclf_clf_synthesis")

fig = plot_synthesis_result(
    quotient_clf,
    verification_clf,
    problem,
    "Induced CLF, verification (∀): " * string(φ),
)
traj, lbl = verification_trajectory(quotient_clf, verification_clf, D_clf)
plot!(fig, ST.Trajectory(traj); label = lbl)
save_fig(fig, "gol_lazar_belta_pclf_clf_verification")

# The induced CLF's slice family under the same regions — the counterpart of the
# `common_slices` figure: the same certified information as the PCLF, folded into one function,
# and its onion is no longer convex. Side by side with their certificate's 11 clean octagons,
# this is the complexity concentration drawn as geometry rather than as a histogram.
fig = plot(; aspect_ratio = :equal, legend = false, size = (700, 620), title = "")
plot!(fig, quotient_clf; what = :slices, show_contours = false, fillalpha = 0.35)
plot!(
    fig,
    problem;
    plot_region = false,
    observation_colors = OBSERVATION_COLORS,
    observation_region_alpha = 0.3,
    observation_linewidth = 2.5,
)
for P in LazySets.array(D_clf)
    plot!(fig, P; fillalpha = 0.0, linecolor = :black, linestyle = :dash, linewidth = 1.5)
end
scatter!(fig, [x0[1]], [x0[2]]; color = :black, markersize = 5)
for (name, R, c) in [
    ("R1", R1, OBSERVATION_COLORS[1]),
    ("R2", R2, OBSERVATION_COLORS[2]),
    ("R3", R3, OBSERVATION_COLORS[3]),
]
    ctr = _centroid(R)
    annotate!(fig, ctr[1], ctr[2], Plots.text(name, 16, c, :center))
end
Dctr_clf = _centroid(LazySets.array(D_clf)[1])
annotate!(fig, Dctr_clf[1], Dctr_clf[2], Plots.text("D", 16, :black, :center))
annotate!(fig, x0[1] - 0.6, x0[2] - 0.7, Plots.text("x₀", 15, :black, :center))
save_fig(fig, "gol_lazar_belta_pclf_clf_slices")

# (1) The problem, poster-style: no domain frame and no legend — each element carries its own
# label — the induced common's largest slice stands in for the working set, and the PCLF's
# synthesized closed loop crosses the scene into its D.
fig = plot(; aspect_ratio = :equal, legend = false, size = (700, 620), title = "")
for slice_list in values(quotient_clf.slices)
    for P in last(slice_list).array
        plot!(fig, P; fillalpha = 0.10, fillcolor = :gray, linealpha = 0.0)
    end
end
plot!(
    fig,
    problem;
    plot_region = false,
    observation_colors = OBSERVATION_COLORS,
    observation_region_alpha = 0.3,
    observation_linewidth = 2.5,
)
for P in LazySets.array(D)
    plot!(fig, P; fillalpha = 0.0, linecolor = :black, linestyle = :dash, linewidth = 1.5)
end
plot!(fig, ST.Trajectory(result.X))
scatter!(fig, [x0[1]], [x0[2]]; color = :black, markersize = 5)
for (name, R, c) in [
    ("R1", R1, OBSERVATION_COLORS[1]),
    ("R2", R2, OBSERVATION_COLORS[2]),
    ("R3", R3, OBSERVATION_COLORS[3]),
]
    ctr = _centroid(R)
    annotate!(fig, ctr[1], ctr[2], Plots.text(name, 16, c, :center))
end
Dctr = _centroid(LazySets.array(D)[1])
annotate!(fig, Dctr[1], Dctr[2], Plots.text("D", 16, :black, :center))
annotate!(fig, x0[1] - 0.6, x0[2] - 0.7, Plots.text("x₀", 15, :black, :center))
save_fig(fig, "gol_lazar_belta_pclf_problem")

# ---------------------------------------------------------
# Augmented 3-D views
# ---------------------------------------------------------
# Everything Plots-based must come before this point — both packages export `plot`.

using CairoMakie

# No PDF for the 3-D views: Cairo writes every mesh triangle as a vector path and the files
# come out above 100 MB. The 3x-resolution PNG is the poster asset.
save_mk(mk, name) = begin
    CairoMakie.save(joinpath(@__DIR__, name * ".png"), mk; px_per_unit = 3)
    return println("wrote ", name, ".png")
end

node_z = Dict((1,) => 0.0, (2,) => 1.0)

# (6) The augmented quotient with one colour per cell — node = the last mode played, so the
# trajectory changes layer exactly when the controller changes mode.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "memory (last mode)",
    zticks = ([0.0, 1.0], ["mode 1", "mode 2"]),
    azimuth = 1.2π,
    elevation = 0.16π,
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
save_mk(mk, "gol_lazar_belta_pclf_pclf_3d_cells")

# (7) Their certificate's by-cell view — one flat layer.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "no memory",
    zticks = ([0.0, 1.0], ["", ""]),
    azimuth = 1.2π,
    elevation = 0.28π,
    title = "Their certificate (one node), one colour per cell",
)
node_z_common = Dict(1 => 0.0)
DI.plot_augmented_bisimulation!(
    ax,
    quotient_common;
    node_z = node_z_common,
    color_by = :state,
    alpha = 0.35,
    show_contours = false,
)
DI.plot_augmented_trajectory!(
    ax,
    quotient_common,
    result_common.X,
    result_common.M;
    node_z = node_z_common,
)
save_mk(mk, "gol_lazar_belta_pclf_common_3d_cells")

# (8) Satisfaction per layer: winning green, losing red.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "memory (last mode)",
    zticks = ([0.0, 1.0], ["mode 1", "mode 2"]),
    azimuth = 1.2π,
    elevation = 0.16π,
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
save_mk(mk, "gol_lazar_belta_pclf_pclf_3d_satisfaction")

# (12) The induced CLF's by-cell view — the same information as the PCLF's, concentrated into
# few fat cells: the complexity redistribution made visible.
mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x1",
    ylabel = "x2",
    zlabel = "no memory",
    zticks = ([0.0, 1.0], ["", ""]),
    azimuth = 1.2π,
    elevation = 0.28π,
    title = "Induced common CLF, one colour per cell",
)
DI.plot_augmented_bisimulation!(
    ax,
    quotient_clf;
    color_by = :state,
    alpha = 0.35,
    show_contours = false,
)
DI.plot_augmented_trajectory!(ax, quotient_clf, result_clf.X, result_clf.M)
save_mk(mk, "gol_lazar_belta_pclf_clf_3d_cells")
