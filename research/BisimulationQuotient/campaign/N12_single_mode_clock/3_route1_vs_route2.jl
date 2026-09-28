# N12 step 2: route 1 against route 2, on the cleanest instance of the geometry-vs-memory trade.
#
# There are TWO baselines here and conflating them is the easy mistake:
#
#   GLB's single node at a FIXED four-facet template.  Certifies 1.134 >= 1 -- no certificate, no
#       sublevel family, nothing to build. That is `2_quotient_on_clock.jl`'s feasibility separation
#       and it is a statement about a template budget, not about route 1.
#   Route 1, the INDUCED common Lyapunov function (`build_common_lyapunov`).  This one works. For the
#       q-cycle over a one-letter alphabet the observer construction collapses to a single subset
#       (reading the only letter maps {0..q-1} onto itself), so
#           V*(x) = max_s ‖R^{-s} x‖_∞,
#       and V*(Ax) = ρ V*(x) exactly. It certifies at the same rate ρ -- with up to 4q facets in one
#       node instead of 4 facets in each of q nodes.
#
# So this experiment is the trade at an IDENTICAL total facet budget: q x 4 distributed against
# 1 x 4q concentrated. Same certificate rate, same guarantee, same problem. Only the geometry differs,
# which is exactly the comparison the paper needs and the one the two-mode benchmarks muddy by also
# changing the rate.
#
# N7 predicts the outcome: `set_difference_decompose` scales as F^2.45 in the plane, so concentrating
# the facets into one node should cost roughly (4q/4)^2.45 = q^2.45 per set difference -- about 50x at
# q = 5 -- against route 2 paying q times as many cheap ones. If route 2 does NOT win here it will not
# win anywhere, because this instance is constructed to isolate the mechanism.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "1_certificate_separation.jl"))

const Q = 5
const RHO = 0.9

f = ST.with_switching(one_mode_system(RHO, Q), HybridSystems.ControlledSwitching())
X = LazySets.Hyperrectangle(; low = [-6.0, -6.0], high = [6.0, 6.0])
R1 = LazySets.Hyperrectangle(; low = [2.0, 0.5], high = [4.0, 2.5])
R2 = LazySets.Hyperrectangle(; low = [-4.0, -2.5], high = [-2.0, -0.5])
problem = PR.BisimulationQuotientProblem(f, X, [R1, R2])

println("x⁺ = ρ R(2π/q) x,  q = ", Q, ", ρ = ", RHO, " -- one mode, no switching\n")

route2 = cycle_pclf(RHO, Q)
route1 = PCLF.build_common_lyapunov(route2)

println(
    "route 2 (path-complete): ",
    length(route2.graph.verts),
    " nodes, rate ",
    route2.JSRapprox,
)
println(
    "route 1 (induced common): ",
    length(route1.graph.verts),
    " node(s), rate ",
    route1.JSRapprox,
)
route1.JSRapprox ≈ route2.JSRapprox || error(
    "the two routes certify different rates ($(route1.JSRapprox) vs $(route2.JSRapprox)); " *
    "`thm:graph-invariance` says they must agree, so this is a defect, not a result",
)
println("same rate, as `thm:graph-invariance` requires -- so the guarantee is identical.\n")

"""
Build one arm and report what it cost, timed warm.

The build is run twice and the faster taken. A quotient construction allocates heavily, so a first
run pays GC and compilation that has nothing to do with the geometry under test -- and this folder has
already retracted one speedup claim (N1) to per-run variance, so a single cold number is not worth
reporting.
"""
function build_arm(name, pclf; atol = 1e-4, max_slices = 30)
    println("--- ", name, " ---")
    best = nothing
    for rep in 1:2
        t0 = time_ns()
        r = build_quotient(
            problem,
            pclf;
            atol = atol,
            level_tol = 1e-2,
            max_slices = max_slices,
            print_level = 0,
        )
        elapsed = (time_ns() - t0) / 1e9
        @printf("  run %d: %.2f s, %d cells\n", rep, elapsed, length(r.quotient.states))
        if best === nothing || elapsed < best.t
            best = (; t = elapsed, quotient = r.quotient)
        end
    end

    q = best.quotient
    parts, faces = PCQ.cell_complexities(q)
    return (;
        name,
        t = best.t,
        cells = length(q.states),
        slices = PCQ.num_slices(q),
        total_faces = sum(faces),
        max_faces = maximum(faces),
        mean_faces = sum(faces) / length(faces),
        total_parts = sum(parts),
        quotient = q,
    )
end

r2 = build_arm("route 2 — $(Q) nodes x 4 facets", route2)
r1 = build_arm("route 1 — 1 node x $(4Q) facets (induced common)", route1)

println("\n", "="^86)
@printf(
    "%-34s %-9s %-9s %-11s %-11s %s\n",
    "arm",
    "cells",
    "slices",
    "Σ faces",
    "max faces",
    "build (s)"
)
println("="^86)
for r in (r2, r1)
    @printf(
        "%-34s %-9d %-9d %-11d %-11d %.2f\n",
        r.name,
        r.cells,
        r.slices,
        r.total_faces,
        r.max_faces,
        r.t
    )
end

@printf(
    "\nroute 2 / route 1:  cells %.2fx   Σ faces %.2fx   build %.2fx\n",
    r2.cells / r1.cells,
    r2.total_faces / r1.total_faces,
    r2.t / r1.t
)
@printf(
    "mean faces per cell: route 2 %.2f, route 1 %.2f  (ratio %.2fx)\n",
    r2.mean_faces,
    r1.mean_faces,
    r1.mean_faces / r2.mean_faces
)

println(
    """

Reading it. A build ratio below 1 means route 2 is faster. The facet budgets are matched by
construction ($(Q) x 4 against 1 x $(4Q)), so any difference is geometry, not certificate strength --
the rates agree to machine precision. N7 measured `set_difference_decompose` at F^2.45 in the plane,
so the per-operation penalty for concentrating facets is the quantity this table prices.""",
)
