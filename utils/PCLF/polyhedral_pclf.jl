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

    G = PCLF.generate_DeBruijn_edges(2, 1; dual = false)

    optimizer = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000)

    θ = π / 6
    R = [
        cos(θ) -sin(θ)
        sin(θ)  cos(θ)
    ]

    # Gmats = :identity
    Gmats = Dict(
        (1,) => [1.0 0.3; 0.0 1.0],   # shear
        (2,) => R,                    # rotation
    )

    pclf = PCLF.compute_polyhedral_pieces_pclf(f, G, optimizer; MLF = true, Gmats = Gmats)
    println("JSR = $(pclf.JSRapprox)")

    gamma = 100.0
    pol1 = PCLF.get_sublevel_set(pclf.pieces[(1,)], gamma)
    pol2 = PCLF.get_sublevel_set(pclf.pieces[(2,)], gamma)
    p = plot(pol1, label="1")
    plot!(pol2, label="2")
    display(p)
end

script()
