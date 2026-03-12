using Dionysos
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PathCompleteFramework

import HybridSystems
import JuMP
import Clarabel
using Plots

function script()
    A1 = [1.5519 0.4474; 7.6412 7.4716]
    A2 = [0.4750 9.1755; 1.8955 0.1850]
    f = HybridSystems.discreteswitchedsystem([A1, A2])
    #G = PCLF.generate_DeBruijn_edges(2, 0, dual=false)
    #println(G)

    G = PCLF.edgeList_to_LabDigraph([
        (1, 3, 2),
        (1, 4, 2),
        (2, 1, 1),
        (2, 2, 1),
        (3, 3, 2),
        (3, 4, 2),
        (4, 1, 1),
        (4, 2, 1),
    ])

    optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

    pclf = PCLF.compute_quadratic_pieces_pclf(f, G, optimizer; MLF = true)
    println("JSR = $(pclf.JSRapprox)")

    gamma = 100.0
    elli1 = PCLF.get_sublevel_set(pclf.pieces[1], gamma)
    elli2 = PCLF.get_sublevel_set(pclf.pieces[2], gamma)
    elli3 = PCLF.get_sublevel_set(pclf.pieces[3], gamma)
    elli4 = PCLF.get_sublevel_set(pclf.pieces[4], gamma)
    p = plot(elli1; label = "1", opacity = 0.5)
    plot!(p, elli2; label = "2", opacity = 0.5)
    plot!(p, elli3; label = "3", opacity = 0.5)
    plot!(p, elli4; label = "4", opacity = 0.5)
    return display(p)
end

script()
