using Dionysos
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PathCompleteFramework

import HybridSystems
using Plots

function script()
    A1 = [1.5519 0.4474; 7.6412 7.4716]
    A2 = [0.4750 9.1755; 1.8955 0.1850]
    f = HybridSystems.discreteswitchedsystem([A1, A2])
    #G = PCLF.generate_DeBruijn_edges(2, 1, dual=false)
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

    pclf = PCLF.compute_quadratic_pieces_pclf(f, G; MLF = true)
    println(pclf.JSRapprox)

    gamma = 100.0
    elli_1 = PCLF.get_sublevel_set(pclf.pieces[1], gamma)
    elli_2 = PCLF.get_sublevel_set(pclf.pieces[2], gamma)
    elli_3 = PCLF.get_sublevel_set(pclf.pieces[3], gamma)
    elli_4 = PCLF.get_sublevel_set(pclf.pieces[4], gamma)
    println(elli_1)
    plot(elli_1)
    plot!(elli_2)
    plot!(elli_3)
    return plot!(elli_4)
end

script()
