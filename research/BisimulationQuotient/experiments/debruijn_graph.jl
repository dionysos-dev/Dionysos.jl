# Bisimulation quotient over a De Bruijn path-complete graph.
#
# Two contracting modes, a De Bruijn graph of order 1 (two nodes), and a rotated template per
# node. Builds the quotient, synthesises a co-safe LTL controller on it, and draws the winning
# region per node and lifted into 3D.

include(joinpath(@__DIR__, "common.jl"))
using Spot

gr()

(; f, problem, R1, R2) = two_mode_problem()

# ---------------------------------------------------------
# Path-complete Lyapunov function
# ---------------------------------------------------------

graph = PCLF.generate_DeBruijn_edges(2, 1; dual = true)
nodes = sort(collect(graph.verts))

pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f,
    graph,
    JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000);
    Gmats = rotation_templates(nodes),
    MLF = true,
    verbose = false,
)
println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

# ---------------------------------------------------------
# Quotient
# ---------------------------------------------------------

(; optimizer, quotient, D) = build_quotient(
    problem,
    pclf;
    atol = 1e-3,
    max_levels = 100,
    max_slices = 8,
    backend = CDDLib.Library(),
)

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

# ---------------------------------------------------------
# Lifted 3D view
# ---------------------------------------------------------
# Loading a Makie backend activates DionysosMakieExt, which provides the lifted recipes.
# GLMakie gives an interactive view; swap in `using CairoMakie` for a static figure.

using GLMakie

node_z = Dict((1,) => 1.0, (2,) => 2.0)

fig = GLMakie.Figure(; size = (900, 700))
ax = GLMakie.Axis3(
    fig[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "node",
    zticks = (collect(values(node_z)), string.(collect(keys(node_z)))),
    title = "Lifted quotient states and closed-loop trajectory",
)
DI.plot_lifted_bisimulation!(
    ax,
    quotient;
    node_z = node_z,
    color_by = :state,
    alpha = 0.2,
    show_contours = false,
)
DI.plot_lifted_trajectory!(ax, quotient, result.X, result.M; node_z = node_z)
display(fig)
