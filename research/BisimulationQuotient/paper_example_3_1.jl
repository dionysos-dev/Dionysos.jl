# Example 3.1 of the paper.
#
# A two-mode system on a square, a common polyhedral PCLF over a one-node graph with a self-loop
# per mode, and a three-part co-safe LTL specification. This is the run the paper's figures and
# reported numbers come from: 10611 quotient states, controllable-set volume 186.728.

include(joinpath(@__DIR__, "common.jl"))
using Spot

gr()

# ---------------------------------------------------------
# System and problem
# ---------------------------------------------------------

A1 = @SMatrix [-0.65 0.32; -0.42 -0.92]
A2 = @SMatrix [0.65 0.32; -0.42 -0.92]
f = HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)])

p = 5.9
X = LazySets.HPolytope([
    LazySets.HalfSpace([1.0, 0.0], p),
    LazySets.HalfSpace([-1.0, 0.0], p),
    LazySets.HalfSpace([0.0, 1.0], p),
    LazySets.HalfSpace([0.0, -1.0], p),
])

R1 = LazySets.HPolytope(
    [-0.9869 -0.1615; -0.0931 0.9957; 0.9659 0.2587; 0.0825 -0.9966],
    [6.6767, 9.2315, 2.3700, -5.9038],
)
R2 = LazySets.HPolytope(
    [0.9993 0.0363; -0.7743 -0.6329; 0.5463 0.8376],
    [-2.1809, 6.3754, -4.8983],
)
R3 = LazySets.HPolytope(
    [-0.9946 -0.1041; 0.5277 0.8494; 0.9999 0.0146; -0.1191 -0.9929],
    [-5.5771, 5.3510, 9.1600, 6.2406],
)

problem = PR.BisimulationQuotientProblem(f, X, [R1, R2, R3])

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

export_optimizer_jld2(optimizer, joinpath(@__DIR__, "paper_example_3_1.jld2"))

fig = plot(; aspect_ratio = :equal)
plot!(fig, quotient; what = :states, node = (1,), show_contours = false)
plot!(fig, problem; opacity = 0.2)
display(fig)

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

fig = plot(; aspect_ratio = :equal, legend = false, title = string(φ))
plot!(
    fig,
    quotient;
    what = :states,
    state_ids = result.controllable_set,
    show_contours = false,
    user_color = :green,
    fillalpha = 1.0,
)
plot!(
    fig,
    quotient;
    what = :states,
    state_ids = result.uncontrollable_set,
    show_contours = false,
    user_color = :red,
    fillalpha = 1.0,
)
plot!(fig, problem; region_alpha = 0.0, observation_region_alpha = 0.0, plot_region = false)
plot!(fig, ST.Trajectory(result.X); label = "Trajectory")
display(fig)
