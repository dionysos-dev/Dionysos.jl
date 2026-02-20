using Dionysos
import HybridSystem
const DI = Dionysos
const UT = DI.Utils
const PCLF = UT.PCLF


function script()
    A1 = [1.5519 0.4474; 7.6412 7.4716]
    A2 = [0.4750 9.1755; 1.8955 0.1850]
    f = HybridSystems.discreteswitchedsystem([A1, A2])
    #G = PCLF.generate_DeBruijn_edges(2, 1, dual=false)
    #println(G)

    G = PCLF.edgeList_to_LabDigraph([(1,3,2), (1,4,2), (2,1,1), (2,2,1), (3,3,2), (3,4,2), (4,1,1), (4,2,1)])

    pclf = PCLF.compute_quadratic_pieces_pclf(f,G, MLF=true)
    print(pclf.JSRapprox)
end