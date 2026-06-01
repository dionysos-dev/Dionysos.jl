using Dionysos
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PathCompleteFramework

import HybridSystems
import JuMP
import Clarabel
using Plots

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

states, trans, alphabet = PCLF.build_observer_graph(graph)
println("Observer states:")
for (k, S) in enumerate(states)
    println("  $k => ", S)
end

println("Transitions:")
for ((k, h), kp) in trans
    println("  ($k, $h) -> $kp")
end

A1 = [1.5519 0.4474; 7.6412 7.4716]
A2 = [0.4750 9.1755; 1.8955 0.1850]
f = HybridSystems.discreteswitchedsystem([A1, A2])

# Quadratic pieces: 
#optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)
#pclf = PCLF.compute_quadratic_pieces_pclf(f, graph, optimizer; MLF = true)

#clf = PCLF.build_common_lyapunov(pclf)
#Vclf = clf.pieces[:clf]

#γ = 1000.0
#xs, ys, vals, mask = PCLF.approximate_sublevel_set(Vclf, γ)

#fig = plot(; aspect_ratio = :equal);

#plot!(PCLF.get_sublevel_set(pclf.pieces[3], γ); label = "3")
#plot!(PCLF.get_sublevel_set(pclf.pieces[1], γ); label = "1")
#plot!(PCLF.get_sublevel_set(pclf.pieces[2], γ); label = "2")
#plot!(PCLF.get_sublevel_set(pclf.pieces[4], γ); label = "4")

#contour!(xs, ys, vals; levels=[γ], aspect_ratio=:equal, linewidth=2)

# Polytopic pieces: 
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
optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)
pclf_poly = PCLF.compute_polyhedral_pieces_pclf(f, graph, optimizer, partitions; MLF = true)

clf_poly = PCLF.build_common_lyapunov(pclf_poly)
Vclf_poly = clf_poly.pieces[:clf]

γ = 0.5
xs_poly, ys_poly, vals_poly, mask_poly = PCLF.approximate_sublevel_set(Vclf_poly, γ)
Sclf = PCLF.approximate_sublevel_set(Vclf_poly, γ)

Sclf_2 = PCLF.get_sublevel_set(Vclf_poly, γ)
#    optimizer_factory = () -> Clarabel.Optimizer()
#)
#println(Sclf_2.parts)

fig = plot(; aspect_ratio = :equal);

#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[3], γ); label = "3")
#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[1], γ); label = "1")
#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[2], γ); label = "2")
#plot!(PCLF.get_sublevel_set(pclf_poly.pieces[4], γ); label = "4")

plot!(Sclf_2; label = "CLF")

#P1 = Sclf_2.parts[1]
#P2 = Sclf_2.parts[2]
#D12 = UT.SemiLinearSet(UT.set_difference_decompose(P2, P1; atol = 1e-6))
#println(D12.parts)

#fig = plot(; aspect_ratio = :equal)
#plot!(P1; alpha = 0.4, label = "P1")
#plot!(P2; alpha = 0.4, label = "P2")
#plot!(D12; alpha = 0.6, label = "P1 \\ P2")
#display(fig)
