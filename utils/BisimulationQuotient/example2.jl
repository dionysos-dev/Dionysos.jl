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

# ---------------------------------------------------------
# Define a stable switched system
# ---------------------------------------------------------
A1 = @SMatrix [
    0.70 0.10;
    0.00 0.65
]

A2 = @SMatrix [
    0.60 -0.15;
    0.10 0.55
]

f = HybridSystems.discreteswitchedsystem([Matrix(A1), Matrix(A2)])

# ---------------------------------------------------------
# Define problem:
# ---------------------------------------------------------

X = Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])
R1 = Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
R2 = Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])
observation_regions = [R1, R2]
problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

# ---------------------------------------------------------
# Compute a common polyhedral PCLF (k = 0 gives 1-node graph)
# ---------------------------------------------------------
G = PCLF.edgeList_to_LabDigraph([
    (1, 1, 1),  # self-loop on node 1 with mode 1
    (1, 2, 2),  # node 1 -> node 2 with mode 2
    (1, 2, 1),  # node 1 -> node 2 with mode 1
    (2, 1, 2),  # node 2 -> node 1 with mode 2
])
sdp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

θ = π / 6
R = [
    cos(θ) -sin(θ)
    sin(θ) cos(θ)
]

θ2 = π / 4
R2 = [
    cos(θ2) -sin(θ2)
    sin(θ2) cos(θ2)
]

# Gmats = :identity
Gmats = Dict(1 => R, 2 => R2)

pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
    f,
    G,
    sdp_optimizer;
    Gmats = Gmats,
    MLF = true,
    verbose = false,
)

println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

# ---------------------------------------------------------
# Instantiate abstraction optimizer
# ---------------------------------------------------------
optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-3)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), 100)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 10)

# ---------------------------------------------------------
# Solve
# ---------------------------------------------------------

MOI.optimize!(optimizer)

construction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
println("Construction time = ", construction_time)

# const FILENAME = joinpath(@__DIR__, "example1.jld2")
# export_optimizer_jld2(optimizer, FILENAME)
# optimizer = import_optimizer_jld2(FILENAME)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))

AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

fig = plot(; layout = (1, 2), aspect_ratio = :equal)

# --- Node 1 ---
plot!(fig[1], bisimulation; what = :states, node = 1, show_contours = false)
title!(fig[1], "Node 1")

# --- Node 2 ---
plot!(fig[2], bisimulation; what = :states, node = 2, show_contours = false)
title!(fig[2], "Node 2")

display(fig)

# ---------------------------------------------------------
# CoSafe LTL control synthesis on the quotient
# ---------------------------------------------------------

using Spot

φ = ltl"F(R1 & F(D))" # ltl"(!R1) U D"
spec = Dionysos.spot_stepper(φ)

x0 = SVector(2.3, 1.5)
_I_ = Hyperrectangle(; low = [x0[1], x0[2]], high = [x0[1], x0[2]])
prob = PR.CoSafeLTLProblem(
    f,
    _I_,
    spec,
    Dict(:D => D, :R1 => R1, :R2 => R2),
    Dict{Symbol, Any}(:D => MP.INNER, :R1 => MP.INNER, :R2 => MP.INNER),
)

optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerCoSafeLTLOnQuotient)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), prob)
MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"), bisimulation)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("ap_to_obs"),
    Dict(:D => -1, :R1 => 1, :R2 => 2),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)
MOI.optimize!(optimizer)

concrete_controller = AB.PCLFBisimulationQuotient.solve_concrete_problem(optimizer)
controllable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("controllable_set"))
uncontrollable_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))

mem0 = AB.PCLFBisimulationQuotient.initial_controller_memory(optimizer, x0)

(X_seq, U_seq, M_seq) = AB.PCLFBisimulationQuotient.simulate_closed_loop(
    f,
    concrete_controller,
    x0,
    mem0;
    N = 50,
)

# ---------------------------------------------------------
# Plot closed-loop trajectory and controlled sets
# ---------------------------------------------------------

fig = plot(; aspect_ratio = :equal);
φ_str = string(φ)
plot!(fig; title = "$φ_str")
plot!(
    bisimulation;
    what = :states,
    state_ids = controllable_set,
    show_contours = false,
    user_color = :green,
    fillalpha = 1.0,
)
plot!(ST.Trajectory(X_seq); label = "Trajectory")
display(fig)

xlims = (-4.5, 4.5)
ylims = (-4.5, 4.5)
fig = plot(; layout = (1, 3), aspect_ratio = :equal, legend = false)
# --- Node 1 ---
title!(fig[1], "Node 1")
xlims!(fig[1], xlims[1], xlims[2])
ylims!(fig[1], ylims[1], ylims[2])
plot!(
    fig[1],
    bisimulation;
    what = :states,
    state_ids = controllable_set,
    node = 1,
    show_contours = false,
    user_color = :green,
    fillalpha = 1.0,
)
plot!(
    fig[1],
    bisimulation;
    what = :states,
    state_ids = uncontrollable_set,
    node = 1,
    show_contours = false,
    user_color = :red,
    fillalpha = 1.0,
)
plot!(fig[1], ST.Trajectory(X_seq); label = "Trajectory")

# --- Node 2 ---
title!(fig[2], "Node 2")
xlims!(fig[2], xlims[1], xlims[2])
ylims!(fig[2], ylims[1], ylims[2])
plot!(
    fig[2],
    bisimulation;
    what = :states,
    state_ids = controllable_set,
    node = 2,
    show_contours = false,
    user_color = :green,
    fillalpha = 1.0,
)
plot!(
    fig[2],
    bisimulation;
    what = :states,
    state_ids = uncontrollable_set,
    node = 2,
    show_contours = false,
    user_color = :red,
    fillalpha = 1.0,
)
plot!(fig[2], ST.Trajectory(X_seq); label = "Trajectory")

# --- All states ---
title!(fig[3], "All states")
xlims!(fig[3], xlims[1], xlims[2])
ylims!(fig[3], ylims[1], ylims[2])
plot!(
    fig[3],
    bisimulation;
    what = :states,
    state_ids = uncontrollable_set,
    show_contours = false,
    user_color = :red,
    fillalpha = 1.0,
)
plot!(
    fig[3],
    bisimulation;
    what = :states,
    state_ids = controllable_set,
    show_contours = false,
    user_color = :green,
    fillalpha = 1.0,
)
plot!(fig[3], ST.Trajectory(X_seq); label = "Trajectory")

display(fig)

# ---------------------------------------------------------
# Plot lifted trajectory and quotient states
# ---------------------------------------------------------

# include("lifted_trajectory_recipes.jl")

# node_z = Dict(1 => 1.0, 2 => 2.0)
# node_colors = Dict(1 => :blue, 2 => :orange)

# fig = GLMakie.Figure(; size = (900, 700))
# ax = GLMakie.Axis3(
#     fig[1, 1];
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

# plot_lifted_trajectory_makie!(ax, bisimulation, X_seq, M_seq; node_z = node_z)

# display(fig)
