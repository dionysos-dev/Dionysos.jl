using StaticArrays
using LinearAlgebra
using JuMP
using Clarabel
using JuMP

import HybridSystems
using LazySets
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

const PCLF = UT.PathCompleteFramework

# Example 3.1 from the paper

# ---------------------------------------------------------
# Define a stable switched system
# ---------------------------------------------------------

A1 = @SMatrix [
    -0.65 0.32;
    -0.42 -0.92
]

A2 = @SMatrix [
    0.65 0.32;
    -0.42 -0.92
]

f = HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)])

# ---------------------------------------------------------
# Define problem:
# ---------------------------------------------------------

p = 5.9
X = Hyperrectangle(; low = [-p, -p], high = [p, p])

H1 = [
    -0.9869 -0.1615
    -0.0931 0.9957
    0.9659 0.2587
    0.0825 -0.9966
]
K1 = [
    6.6767
    9.2315
    2.3700
    -5.9038
]
R1 = HPolytope(H1, K1)

H2 = [
    0.9993 0.0363
    -0.7743 -0.6329
    0.5463 0.8376
]
K2 = [
    -2.1809
    6.3754
    -4.8983
]
R2 = HPolytope(H2, K2)

H3 = [
    -0.9946 -0.1041
    0.5277 0.8494
    0.9999 0.0146
    -0.1191 -0.9929
]
K3 = [
    -5.5771
    5.3510
    9.1600
    6.2406
]
R3 = HPolytope(H3, K3)

observation_regions = [R1, R2, R3]
problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

# ---------------------------------------------------------
# Common polyhedral PCLF
# ---------------------------------------------------------
#L = [
#    -0.0625 1.0
#     0.6815 1.0
#     0.9947 0.6868
#     0.9947 -0.0678
# ]

#rho = 0.94
M = 2   # two modes

# Common Lyapunov graph = one node with one self-loop per mode
graph = PCLF.generate_DeBruijn_edges(M, 1)

G = PCLF.edgeList_to_LabDigraph([(1, 1, 1), (2, 1, 1), (2, 1, 2), (1, 2, 2)])

# Recover the unique vertex
# v = first(graph.verts)
# println("Vertex: ", v)
# One piece: V(x) = ||Lx||_inf
# piece = PCLF.PolyhedralPiece(L, ones(size(L, 1)))

# PCLF with one piece
# pclf = PCLF.PCLF(graph, Dict(v => piece), rho)

v1 = [1.0, 0.0]
v2 = [1.0, 1.0]
v3 = [0.0, 1.0]
v4 = [-1.0, 1.0]
v5 = [-1.0, 0.0]

C1 = hcat(v1, v2)
C2 = hcat(v2, v3)
C3 = hcat(v3, v4)
C4 = hcat(v4, v5)

partitions = Dict((1,) => [C1, C2, C3, C4], (2,) => [C1, C2, C3, C4])

optimizer_pclf = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)
pclf = PCLF.compute_polyhedral_pieces_pclf(f, graph, optimizer_pclf, partitions; MLF = true)

println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

# ---------------------------------------------------------
# Instantiate abstraction optimizer
# ---------------------------------------------------------

optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf)
MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-4)
MOI.set(optimizer, MOI.RawOptimizerAttribute("level_tol"), 1e-2)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), 100)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("ΓX"), 10.0) # nothing
# MOI.set(optimizer, MOI.RawOptimizerAttribute("nb_levels"), 12) # nothing
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 50)

# ---------------------------------------------------------
# Solve
# ---------------------------------------------------------
MOI.optimize!(optimizer)
construction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
println("Construction time = ", construction_time)

# const FILENAME = joinpath(@__DIR__, "example_3_1.jld2")
#AB.PCLFBisimulationQuotient.export_optimizer_jld2(optimizer, FILENAME)
# optimizer = AB.PCLFBisimulationQuotient.import_optimizer_jld2(FILENAME)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))

AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

fig = plot(; aspect_ratio = :equal);
# fig = plot();
# plot!(bisimulation; what = :slices, show_contours = false)
plot!(bisimulation; what = :states, by = :states, node = (2,), show_contours = false)
plot!(problem; opacity = 0.2)
display(fig)

# ---------------------------------------------------------
# CoSafe LTL control synthesis on the quotient
# ---------------------------------------------------------

#using Spot

#φ = ltl"((!R2 U D) & F(R1) & ((R3 -> X(!R1)) U D))"

#x0 = SVector(-6.0, 7.5)
#x0 = SVector(9.0, 0.0)     # initial point I in the paper
#x0 = SVector(-4.0, -7.0)   # initial point a in the paper
#_I_ = Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]])

#prob = PR.CoSafeLTLProblem(
#    f,
#    _I_,
#    φ,
#    Dict(:D => D, :R1 => R1, :R2 => R2, :R3 => R3), # no really useful since we have ap_to_obs, but let's be explicit
#    Dict{Symbol, Any}(:D => MP.INNER, :R1 => MP.INNER, :R2 => MP.INNER, :R3 => MP.INNER), # no really useful, but let's be explicit
#    true,
#)

#optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), prob)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"), bisimulation)
#MOI.set(
#    optimizer,
#    MOI.RawOptimizerAttribute("ap_to_obs"),
#    Dict(:D => -1, :R1 => 1, :R2 => 2, :R3 => 3),
#)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)
#MOI.optimize!(optimizer)

#concrete_controller = AB.PCLFBisimulationQuotient.solve_concrete_problem(optimizer)
#controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
#uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))

#Vctrl = AB.PCLFBisimulationQuotient.get_volume(bisimulation, controllable_set)
#println("Volume of controllable set = ", Vctrl)

#mem0 = AB.PCLFBisimulationQuotient.initial_controller_memory(optimizer, x0)

#(X_seq, U_seq, M_seq) = AB.PCLFBisimulationQuotient.simulate_closed_loop(
#    f,
#    concrete_controller,
#    x0,
#    mem0;
#    N = 50,
#)

# ---------------------------------------------------------
# Plot closed-loop trajectory and controlled sets
# ---------------------------------------------------------

#fig = plot(; aspect_ratio = :equal, legend = false)
#φ_str = string(φ)
#plot!(fig; title = "$φ_str")
#plot!(
#    bisimulation;
#    what = :states,
#    state_ids = controllable_set,
#    show_contours = false,
#    user_color = :green,
#    fillalpha = 1.0,
#)
# plot!(bisimulation; what = :states, state_ids = uncontrollable_set, show_contours = false, user_color = :red, fillalpha = 1.0)
#plot!(problem; region_alpha = 0.0, observation_region_alpha = 0.0, plot_region = false)
#plot!(ST.Trajectory(X_seq); label = "Trajectory")
#display(fig)

# ---------------------------------------------------------
# Plot lifted trajectory and quotient states
# ---------------------------------------------------------

# include("lifted_trajectory_recipes.jl")

# node_z = Dict(1 => 1.0)
# node_colors = Dict(1 => :blue)

# fig = GLMakie.Figure(size = (900, 700))
# ax = GLMakie.Axis3(
#     fig[1, 1],
#     xlabel = "x₁",
#     ylabel = "x₂",
#     zlabel = "node",
#     zticks = (collect(values(node_z)), string.(collect(keys(node_z)))),
#     title = "Lifted quotient states and closed-loop trajectory",
# )

# plot_lifted_bisimulation_makie!(
#     ax,
#     bisimulation;
#     node_z = node_z,
#     color_by = :state,
#     alpha = 0.2,
#     show_contours = false,
# )

# plot_lifted_trajectory_makie!(
#     ax,
#     bisimulation,
#     X_seq,
#     M_seq;
#     node_z = node_z,
# )

# display(fig)
