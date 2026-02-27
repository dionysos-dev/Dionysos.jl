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

    G = PCLF.generate_DeBruijn_edges(2, 2; dual = false)

    optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

    pclf = PCLF.compute_polyhedral_pieces_pclf(f, G, optimizer; MLF = true)
    return println(pclf.JSRapprox)
end

script()
