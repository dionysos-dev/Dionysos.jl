using Dionysos
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PathCompleteFramework

import HybridSystems
import JuMP
import Clarabel
using Plots

function script()
    #A1 = @SMatrix [
    #    -0.65 0.32;
    #    -0.42 -0.92
    #]

    #A2 = @SMatrix [
    #    0.65 0.32;
    #    -0.42 -0.92   
    #]

    f = HybridSystems.discreteswitchedsystem([A1, A2])

    G = PCLF.generate_DeBruijn_edges(2, 2; dual = false)

    #G = PCLF.edgeList_to_LabDigraph([
    #    (1, 1, 1),
    #    (1, 2, 1),
    #    (1, 2, 2),
    #    (2, 1, 2)
    #])

    optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

    θ = π / 6
    R = [
        cos(θ) -sin(θ)
        sin(θ) cos(θ)
    ]

    Gmats = :identity
    #Gmats = Dict(
    #    (1,) => [1.0 0.3; 0.0 1.0],   # shear
    #    (2,) => R,                    # rotation
    #)
    #Gmats = [[-0.6672110938517144 0.4951760688812728; 0.7397289852833637 0.36222336963570656], [-1.115084393756029 0.6782490219162379; -0.7422010232488797 -0.6925104522392596]]
    #Gmats = [randn(2, 2) for _ in 1:4]
    #println(Gmats)

    pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(f, G, optimizer; MLF = true, Gmats = Gmats)
    println("JSR = $(pclf.JSRapprox)")

    #gamma = 100.0
    #pol1 = PCLF.get_sublevel_set(pclf.pieces[(1,)], gamma)
    #pol2 = PCLF.get_sublevel_set(pclf.pieces[(2,)], gamma)
    #p = plot(pol1; label = "1")
    #plot!(pol2; label = "2")
    #return display(p)
end

function new_polyhedral()

    A1 = [-0.65 0.32; -0.42 -0.92]
    A2 = [0.65 0.32; -0.42 -0.92]
    f = HybridSystems.discreteswitchedsystem([A1, A2])

    G = PCLF.generate_DeBruijn_edges(2, 1; dual = false)

    #G = PCLF.edgeList_to_LabDigraph([
    #    (1, 1, 1),
    #    (2, 1, 1),
    #    (2, 1, 2),
    #    (1, 2, 2)
    #])
    println(G)

    # standard basis
    #e1 = [1.0, 0.0]
    #e2 = [0.0, 1.0]
    # define cones (each column = one ray)
    #C1_1 = hcat( e1,  e2)   # Q1
    #C2_1 = hcat(-e1,  e2)   # Q2
    #C3_1 = hcat(-e1, -e2)   # Q3
    #C4_1 = hcat( e1, -e2)   # Q4

    #a = 1 / sqrt(2)
    #r1 = [ a,  a]   # 45 degrees
    #r2 = [-a,  a]   # 135 degrees
    #r3 = [-a, -a]   # 225 degrees
    #r4 = [ a, -a]   # 315 degrees

    #C1_2 = hcat(r1, r2)
    #C2_2 = hcat(r2, r3)
    #C3_2 = hcat(r3, r4)
    #C4_2 = hcat(r4, r1)

    v1 = [ 1.0,  0.0]
    v2 = [ 1.0,  1.0]
    v3 = [ 0.0,  1.0]
    v4 = [-1.0,  1.0]
    v5 = [-1.0,  0.0]

    C1 = hcat(v1, v2)
    C2 = hcat(v2, v3)
    C3 = hcat(v3, v4)
    C4 = hcat(v4, v5)

    partitions = Dict(
        (1,) => [C1, C2, C3, C4],
        (2,) => [C1, C2, C3, C4],
    )

    optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

    pclf = PCLF.compute_polyhedral_pieces_pclf(f, G, optimizer, partitions; MLF = true)
    println("JSR = $(pclf.JSRapprox)")

    gamma = 10.0
    pol1 = PCLF.get_sublevel_set(pclf.pieces[(1,)], gamma)
    pol2 = PCLF.get_sublevel_set(pclf.pieces[(2,)], gamma)
    p = plot(pol1; label = "1")
    plot!(pol2; label = "2")
    return display(p)

end 

#script()
new_polyhedral()
