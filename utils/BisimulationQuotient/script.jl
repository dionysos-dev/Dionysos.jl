using StaticArrays
using LinearAlgebra
using JuMP
using Clarabel
using MathOptInterface

import HybridSystems
using LazySets
using Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
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
#  Define region X, terminal region D, and observation regions
# ---------------------------------------------------------

# X = [-2,2]^2
X = Hyperrectangle(; low = [-2.0, -2.0], high = [2.0, 2.0])

# D = [-0.2,0.2]^2
D = Hyperrectangle(; low = [-0.2, -0.2], high = [0.2, 0.2])

# Two observation regions
R1 = Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
R2 = Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])

observation_regions = [R1, R2]

problem = PR.BisimulationQuotientProblem(f, X, D, observation_regions)

# ---------------------------------------------------------
# Compute a common polyhedral PCLF (k = 0 gives 1-node graph)
# ---------------------------------------------------------
G = PCLF.generate_DeBruijn_edges(2, 0; dual = false)

sdp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

pclf = PCLF.compute_polyhedral_pieces_pclf(f, G, sdp_optimizer; MLF = true, verbose = false)

println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

# ---------------------------------------------------------
# Define Lyapunov levels Γ
# Must be increasing, with Γ[i] ≈ γ * Γ[i+1] structure
# ---------------------------------------------------------
# γ = pclf.JSRapprox
# Γ0 = 200.0
# N = 8
# Γ = [Γ0 / (γ^i) for i in 0:(N-1)]

# println("Levels Γ = ", Γ)

# ---------------------------------------------------------
# Instantiate abstraction optimizer
# ---------------------------------------------------------
optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf)
MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-2)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("Γ"), Float64.(Γ))
MOI.set(optimizer, MOI.RawOptimizerAttribute("num_levels"), 5)

# ---------------------------------------------------------
# Solve
# ---------------------------------------------------------
MOI.optimize!(optimizer)

println(
    "Construction time = ",
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec")),
)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
obs_partition = MOI.get(optimizer, MOI.RawOptimizerAttribute("obs_partition"))

println("Number of abstract states = ", length(bisimulation.states))

fig = plot(; aspect_ratio = :equal);
# plot!(problem)                # problem geometry
# plot!(obs_partition)          # observation partition
plot!(bisimulation; what = :slices, mode = 1)
# plot!(bisimulation; what = :states, mode = 1, by = :obs)
# plot!(bisimulation; what = :states, mode = 1, by = :slice)
display(fig)
