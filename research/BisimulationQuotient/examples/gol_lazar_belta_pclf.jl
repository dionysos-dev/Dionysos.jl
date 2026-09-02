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

include(joinpath(@__DIR__, "common.jl"))
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

display(plot_quotient(quotient, problem, "Quotient cells of node (1,)"; node = (1,)))

# ---------------------------------------------------------
# Co-safe LTL synthesis
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

# ---------------------------------------------------------
# Closed-loop trajectory over the winning region
# ---------------------------------------------------------

fig = plot_synthesis_result(quotient, result, problem, string(φ))
display(fig)
