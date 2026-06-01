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
p = 1.7
X = Hyperrectangle(; low = [-p, -p], high = [p, p])
R1 = Hyperrectangle(; low = [0.8, 0.8], high = [1.5, 1.5])
#R2 = Hyperrectangle(; low = [-1.5, 0.8], high = [-0.8, 1.5])
observation_regions = [R1]
problem = PR.BisimulationQuotientProblem(f, X, observation_regions)

# -- Computation of a PCLF: 

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

fig = plot(; aspect_ratio = :equal);
plot!(PCLF.get_sublevel_set(pclf_poly.pieces[3], 0.5); label = "3")
plot!(PCLF.get_sublevel_set(pclf_poly.pieces[1], 0.5); label = "1")
plot!(PCLF.get_sublevel_set(pclf_poly.pieces[2], 0.5); label = "2")
plot!(PCLF.get_sublevel_set(pclf_poly.pieces[4], 0.5); label = "4")
display(fig)

# -- Computation of the abtraction: 

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

# ---------------------------------------------------------
# Solve
# ---------------------------------------------------------
MOI.optimize!(optimizer)
construction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
println("Construction time = ", construction_time)

bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))

AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

fig = plot(; aspect_ratio = :equal);
# fig = plot();
# plot!(bisimulation; what = :slices, show_contours = false)
plot!(bisimulation; what = :states, by = :states, node = 2, show_contours = false)
plot!(problem; opacity = 0.2)
display(fig)

# -- The induced CLF -----------------------------------------------------------------------------------------

#states, trans, alphabet = PCLF.build_observer_graph(graph)
#clf_poly = PCLF.build_common_lyapunov(pclf_poly)

#Vclf_poly = clf_poly.pieces[:clf]

#Sclf_2 = PCLF.get_sublevel_set(Vclf_poly, 0.5)

#fig = plot(; aspect_ratio = :equal);

#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[3], 0.5); label = "3")
#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[1], 0.5); label = "1")
#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[2], 0.5); label = "2")
#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[4], 0.5); label = "4")

#plot!(Sclf_2; label = "CLF")

# Computation of the abstraction

#levels, D = AB.PCLFBisimulationQuotient.build_levels_and_terminal_set(
#    clf_poly,   # <-- CLF, not pclf_poly
#    X,
#    observation_regions;
#    tol = 1e-3,
#    max_levels = 50,
#)

#fig = plot(; aspect_ratio = :equal);
#for τ in levels
#    Sτ = PCLF.get_sublevel_set(Vclf_poly, τ)
#    plot!(fig, Sτ; alpha = 0.05, label = nothing)
#end
#plot!(X)
#plot!(D)
#display(fig)

#optimizer = MOI.instantiate(AB.PCLFBisimulationQuotient.OptimizerBisimulationQuotient)

#MOI.set(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient_problem"), problem)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("pclf"), clf_poly)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), true)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("atol"), 1e-4)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("level_tol"), 1e-2)
#MOI.set(optimizer, MOI.RawOptimizerAttribute("max_slices"), 20)

#MOI.optimize!(optimizer)
#construction_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("construction_time_sec"))
#println("Construction time = ", construction_time)

#bisimulation = MOI.get(optimizer, MOI.RawOptimizerAttribute("bisimulation_quotient"))
#D = MOI.get(optimizer, MOI.RawOptimizerAttribute("D"))

#AB.PCLFBisimulationQuotient.print_bisimulation_stats(bisimulation)

#fig = plot(; aspect_ratio = :equal);
# fig = plot();
# plot!(bisimulation; what = :slices, show_contours = false)
#plot!(bisimulation; what = :states, by = :states, node = :clf, show_contours = false)
#plot!(problem; opacity = 0.2)
#display(fig)
