# A worked example for the paper: the lifted quotient and a closed-loop trajectory that uses the
# memory.
#
# The point of the figure is to make the lifting visible. On the cycle graph the node advances or
# retreats with every mode, so a satisfying run threads through all the layers rather than staying
# on one — which is what distinguishes the lifted abstraction from a quotient of the state space
# alone.
#
#     A₁ = ρ R(2π/q),   A₂ = ρ R(-2π/q)
#     graph  ℤ_q  with edges (s, s+1, mode 1) and (s, s-1, mode 2)
#     PCLF   V_s(x) = ‖R^{-s}x‖_∞     -- written in closed form, rate exactly ρ
#
# ρ is chosen close to 1 so that trajectories take many steps before reaching the terminal set; a
# fast contraction would reach it in three or four steps and show almost nothing.

include(joinpath(@__DIR__, "common.jl"))

using Spot
using Printf

const Q = 5              # number of graph nodes, hence layers in the lifted figure
const RHO = 0.94         # slow contraction, so the trajectory is long enough to be legible

rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

rotation_system(ρ, q) =
    HybridSystems.discreteswitchedsystem([ρ .* rotation(2π / q), ρ .* rotation(-2π / q)])

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

graph, pclf = cycle_pclf(RHO, Q)
f = rotation_system(RHO, Q)

println(
    "graph on ℤ_$Q: complete = ",
    PCLF.is_complete(graph, 1:2),
    ", deterministic = ",
    PCLF.is_deterministic(graph, 1:2),
)
println("certificate rate = ", pclf.JSRapprox)

X = LazySets.Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
R1 = LazySets.Hyperrectangle(; low = [0.9, 0.7], high = [1.7, 1.5])
R2 = LazySets.Hyperrectangle(; low = [-1.7, 0.7], high = [-0.9, 1.5])
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])

(; quotient, D) =
    build_quotient(problem, pclf; atol = 1e-3, max_slices = 9, print_level = 0)

# ---------------------------------------------------------
# Synthesis
# ---------------------------------------------------------
# Reaching R1 before the terminal set takes the trajectory around the origin. On the cycle graph
# each mode advances or retreats the node, so the run threads through all five layers: from this
# initial state it takes ten steps and changes node nine times.

φ = ltl"F(R1 & F(D))"
x0 = SVector(-1.0, 1.7)

result = synthesize_cosafe_ltl(
    f,
    quotient,
    Dionysos.spot_stepper(φ),
    Dict(:D => D, :R1 => R1, :R2 => R2),
    Dict(:D => -1, :R1 => 1, :R2 => 2),
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
# Planar view
# ---------------------------------------------------------

fig = plot(; aspect_ratio = :equal, legend = false, title = string(φ))
plot!(
    fig,
    quotient;
    what = :states,
    state_ids = result.controllable_set,
    show_contours = false,
    user_color = :seagreen,
    fillalpha = 0.9,
)
plot!(
    fig,
    quotient;
    what = :states,
    state_ids = result.uncontrollable_set,
    show_contours = false,
    user_color = :indianred,
    fillalpha = 0.9,
)
plot!(
    fig,
    problem;
    region_alpha = 0.0,
    observation_region_alpha = 0.35,
    plot_region = false,
)
plot!(fig, ST.Trajectory(result.X); label = "trajectory")
savefig(fig, joinpath(@__DIR__, "lifted_showcase_planar.png"))
println("wrote lifted_showcase_planar.png")

# ---------------------------------------------------------
# Lifted view
# ---------------------------------------------------------
# Loading a Makie backend activates DionysosMakieExt. CairoMakie gives a static figure suitable
# for the paper; GLMakie gives an interactive one.

using CairoMakie

node_z = Dict(s => Float64(s) for s in 0:(Q - 1))

mk = CairoMakie.Figure(; size = (900, 700))
ax = CairoMakie.Axis3(
    mk[1, 1];
    xlabel = "x₁",
    ylabel = "x₂",
    zlabel = "graph node",
    zticks = (collect(values(node_z)), string.(collect(keys(node_z)))),
    title = "Lifted quotient and closed-loop trajectory",
)

# Every cell is drawn, not only the winning ones, so each layer shows the whole partition carried
# by that node and the lifting is legible as five stacked copies of the state space.
DI.plot_lifted_bisimulation!(
    ax,
    quotient;
    node_z = node_z,
    color_by = :node,
    node_colors = Dict(
        s => c for
        (s, c) in zip(0:(Q - 1), [:steelblue, :seagreen, :goldenrod, :orchid, :tomato])
    ),
    alpha = 0.20,
    show_contours = false,
)
DI.plot_lifted_trajectory!(ax, quotient, result.X, result.M; node_z = node_z)

CairoMakie.save(joinpath(@__DIR__, "lifted_showcase_3d.png"), mk)
println("wrote lifted_showcase_3d.png")
