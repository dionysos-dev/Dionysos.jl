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

X = Hyperrectangle(; low = [-10.0, -10.0], high = [10.0, 10.0])
R1 = Hyperrectangle(; low = [0.5, 0.5], high = [0.9, 0.9])
R2 = Hyperrectangle(; low = [0.5, -0.5], high = [0.9, -0.1])
observation_regions = [R1, R2]
problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

# ---------------------------------------------------------
# Compute a common polyhedral PCLF (k = 0 gives 1-node graph)
# ---------------------------------------------------------

G = PCLF.generate_DeBruijn_edges(2, 2; dual = false)
sdp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

θ = π / 6
R = [
    cos(θ) -sin(θ)
    sin(θ) cos(θ)
]

Gmats = :identity
# Gmats = Dict(
#     (1,1) => R,   # rotation
#     (2,1) => R,   # rotation
#     (1,2) => R,   # rotation
#     (2,2) => R,   # rotation    
# )

pclf = PCLF.compute_polyhedral_pieces_pclf(
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
MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-3)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), 100)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 8)

# ---------------------------------------------------------
# Solve
# ---------------------------------------------------------
MOI.optimize!(optimizer)
construction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
println("Construction time = ", construction_time)

# const FILENAME = joinpath(@__DIR__, "example2.jld2")
# AB.PCLFBisimulationQuotient.export_optimizer_jld2(optimizer, FILENAME)
# optimizer = AB.PCLFBisimulationQuotient.import_optimizer_jld2(FILENAME)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))

AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

fig = plot(; aspect_ratio = :equal);
plot!(problem)
plot!(bisimulation; what = :states, node = (1,), by = :obs)
display(fig)
