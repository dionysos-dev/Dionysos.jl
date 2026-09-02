# Bisimulation quotient over a hand-built path-complete graph.
#
# Same two-mode system as `debruijn_graph.jl`, but the path-complete graph is written out edge
# by edge instead of generated, and the two nodes carry templates rotated by different angles.
# The point of the comparison is that a graph chosen by hand can be smaller than the De Bruijn
# graph of the same order and still be path-complete.

include(joinpath(dirname(@__DIR__), "common.jl"))
using Spot

gr()

(; f, problem, R1, R2) = two_mode_problem()

# ---------------------------------------------------------
# Path-complete Lyapunov function
# ---------------------------------------------------------

graph = PCLF.edgeList_to_LabDigraph([
    (1, 1, 1),  # self-loop on node 1 with mode 1
    (1, 2, 2),  # node 1 -> node 2 with mode 2
    (1, 2, 1),  # node 1 -> node 2 with mode 1
    (2, 1, 2),  # node 2 -> node 1 with mode 2
])

# A distinct rotation per node, so this is not the uniform template of `debruijn_graph.jl`.
rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
templates = Dict(1 => rotation(π / 6), 2 => rotation(π / 4))

pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f,
    graph,
    JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000);
    Gmats = templates,
    MLF = true,
    verbose = false,
)
println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

# ---------------------------------------------------------
# Quotient
# ---------------------------------------------------------

(; optimizer, quotient, D) =
    build_quotient(problem, pclf; atol = 1e-3, max_levels = 100, max_slices = 10)

display(plot_cells_by_node(quotient))

# ---------------------------------------------------------
# Co-safe LTL synthesis
# ---------------------------------------------------------

φ = ltl"F(R1 & F(D))"
x0 = SVector(2.3, 1.5)

result = synthesize_cosafe_ltl(
    f,
    quotient,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2),
    Dict(:D => -1, :R1 => 1, :R2 => 2),
    x0,
)

# ---------------------------------------------------------
# Figures
# ---------------------------------------------------------

display(plot_controllable_set(quotient, result.controllable_set, result.X, string(φ)))
display(
    plot_partition_by_node(
        quotient,
        result.controllable_set,
        result.uncontrollable_set,
        result.X,
    ),
)
