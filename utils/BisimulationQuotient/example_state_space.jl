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

A1 = (1.0/10.0)*[1.5519 0.4474; 7.6412 7.4716]
A2 = (1.0/10.0)*[0.4750 9.1755; 1.8955 0.1850]
f = HybridSystems.discreteswitchedsystem([A1, A2])

# -- Definition of the problem:
p = 2.5

θ = deg2rad(10.0)

R = [
    cos(θ) -sin(θ);
    sin(θ) cos(θ)
]

normals = [R * [1.0, 0.0], R * [-1.0, 0.0], R * [0.0, 1.0], R * [0.0, -1.0]]

# State space
#= X = HPolytope([
    HalfSpace(normals[1], p),
    HalfSpace(normals[2], p),
    HalfSpace(normals[3], p),
    HalfSpace(normals[4], p),
]) =#
X = HPolytope([
    HalfSpace([1.0, 0.0], 3.1),   # x ≤ 3.1
    HalfSpace([-1.0, 0.0], 3.2),   # x ≥ -3.2
    HalfSpace([0.0, 1.0], 3.6),   # y ≤ 3.6
    HalfSpace([0.0, -1.0], 2.9),   # y ≥ -2.9
    HalfSpace([0.8, 0.6], 2.7),
    HalfSpace([-0.8, -0.6], 2.7),
    HalfSpace([0.6, -0.8], 2.8),
    HalfSpace([-0.6, 0.8], 2.8),
])

# R1 = [0.8, 1.5] × [0.8, 1.5]
R1 = HPolytope([
    HalfSpace([1.0, 0.0], 1.5),
    HalfSpace([-1.0, 0.0], -0.8),
    HalfSpace([0.0, 1.0], 1.5),
    HalfSpace([0.0, -1.0], -0.8),
])

# R2 = [-1.5, -0.8] × [0.8, 1.5]
R2 = HPolytope([
    HalfSpace([1.0, 0.0], -0.8),
    HalfSpace([-1.0, 0.0], 1.5),
    HalfSpace([0.0, 1.0], 1.5),
    HalfSpace([0.0, -1.0], -0.8),
])

observation_regions = [R1]
problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

graph = PCLF.edgeList_to_LabDigraph([
    (1, 2, 1),
    (2, 1, 1),
    (2, 4, 1),
    (2, 3, 1),
    (3, 4, 1),
    (4, 3, 2),
    (4, 4, 2),
    (4, 1, 2),
])

v1 = [1.0, 0.0]
v2 = [1.0, 1.0]
v3 = [0.0, 1.0]
v4 = [-1.0, 1.0]
v5 = [-1.0, 0.0]

C1 = hcat(v1, v2)
C2 = hcat(v2, v3)
C3 = hcat(v3, v4)
C4 = hcat(v4, v5)

partitions = Dict(
    1 => [C1, C2, C3, C4],
    2 => [C1, C2, C3, C4],
    3 => [C1, C2, C3, C4],
    4 => [C1, C2, C3, C4],
)
lp_optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)
pclf_poly =
    PCLF.compute_polyhedral_pieces_pclf(f, graph, lp_optimizer, partitions; MLF = true)
println("Computed JSR upper bound / contraction rate = ", pclf_poly.JSRapprox)

optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), pclf_poly)
MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-4)
MOI.set(optimizer, MOI.RawOptimizerAttribute("level_tol"), 1e-2)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("max_levels"), 100)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("ΓX"), 10.0) # nothing
# MOI.set(optimizer, MOI.RawOptimizerAttribute("nb_levels"), 12) # nothing
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 30)

MOI.optimize!(optimizer)
construction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
println("Construction time = ", construction_time)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))

AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

fig = plot(; aspect_ratio = :equal);
plot!(bisimulation; what = :states, by = :states, node = 1, show_contours = false)
plot!(problem; opacity = 0.2)
display(fig)
