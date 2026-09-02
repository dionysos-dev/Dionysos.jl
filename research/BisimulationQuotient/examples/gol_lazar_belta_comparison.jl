# Our PCLF machinery against the hand-crafted certificate of Gol, Ding, Lazar & Belta, on their
# own Example 3.1 — at their own budget.
#
# Their certificate is one node carrying four template rows (`V = ‖Lx‖_∞`, eight half-space
# facets). The fair question is what our solver does with the *same* per-node budget: the conic
# solver at partition order 1 — four cones, hence four rows per node — with the partition's
# orientation searched over a grid, once without memory (|S| = 1) and once with (the order-1
# De Bruijn graph, |S| = 2). Both arms get the same orientation search — a single angle shared by
# all nodes — because giving one arm a richer parameterisation manufactured false separations
# twice during this project. (A first draft used the symmetric 2n-face template and certified
# nothing at any angle: that family has two rows, half their budget — the mismatch was the
# result.)
#
# What is comparable and what is not, per the folder README: the certified contraction rate is a
# bound on the joint spectral radius and compares across certificates; raw cell counts do not,
# because each certificate tiles its own region — so cells are reported *with* the covered volume
# and per covered volume. Rates are certified exactly (`certify_pclf`, linear programming), not
# read off the solver's bisection.

include(joinpath(dirname(@__DIR__), "common.jl"))

using Printf
import HiGHS

(; f, problem, X, R1, R2, R3) = gol_lazar_belta_problem()
const L = gol_lazar_belta_L()

exact_rate(pclf) = PCLF.certify_pclf(
    PCLF.PCLF(pclf.graph, pclf.pieces, pclf.JSRapprox),
    f,
    HiGHS.Optimizer,
).rate

# ---------------------------------------------------------
# The three certificates
# ---------------------------------------------------------

pclf_theirs = PCLF.PCLF(
    PCLF.generate_DeBruijn_edges(2, 0),
    Dict{Any, PCLF.AbstractPiece}(1 => PCLF.PolyhedralPiece(L, ones(size(L, 1)))),
    0.94,
)

rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

# One angle, shared by every node, rotating the order-1 conic partition — four cones, four rows
# per node, their budget. The partition is symmetric under a quarter turn, so [0, π/2) covers the
# orientations once. Only a certified contraction (exact rate < 1) counts as found.
function best_conic(k; nθ = 24)
    graph = PCLF.generate_DeBruijn_edges(2, k)
    sdp = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 1000,
        MOI.Silent() => true,
    )
    base = PCLF.conic_partitions_2d(1)
    best = nothing
    best_rate = Inf
    best_θ = NaN
    for θ in range(0, π / 2; length = nθ + 1)[1:(end - 1)]
        cones = [rotation(θ) * C for C in base]
        partitions = Dict(v => cones for v in graph.verts)
        pclf = PCLF.compute_polyhedral_pieces_pclf(f, graph, sdp, partitions; MLF = true)
        isfinite(pclf.JSRapprox) || continue
        rate = exact_rate(pclf)
        if rate < 1.0 && rate < best_rate
            best, best_rate, best_θ = pclf, rate, θ
        end
    end
    best === nothing &&
        error("No contracting certificate found for k = $k at any orientation.")
    # The quotient's slice recursion runs on `JSRapprox`; give it the certified rate, not the
    # bisection's upper end.
    return PCLF.PCLF(best.graph, best.pieces, best_rate), best_rate, best_θ
end

t = @elapsed begin
    global pclf_k0, rate_k0, θ_k0 = best_conic(0)
    global pclf_k1, rate_k1, θ_k1 = best_conic(1)
end
rate_theirs = exact_rate(pclf_theirs)

println("certificate search: $(round(t; digits = 1)) s for both arms")
@printf("their hand-crafted L       : ρ = %.6f  (1 node × 4 facets)\n", rate_theirs)
@printf("our solver, no memory (k=0): ρ = %.6f  at θ = %.3f\n", rate_k0, θ_k0)
@printf("our solver, memory    (k=1): ρ = %.6f  at θ = %.3f\n", rate_k1, θ_k1)

# ---------------------------------------------------------
# The three quotients, identically parameterised
# ---------------------------------------------------------

function quotient_stats(name, pclf)
    t = @elapsed (; quotient, D) =
        build_quotient(problem, pclf; atol = 1e-3, max_slices = 40, print_level = 0)
    parts, faces = PCQ.cell_complexities(quotient)
    covered = PCQ.get_volume(quotient, keys(quotient.states); backend = CDDLib.Library())
    nslices = length(first(values(quotient.slices)))
    return (;
        name,
        quotient,
        D,
        t,
        ncells = length(quotient.states),
        nslices,
        Σfacets = sum(faces),
        maxfacets = maximum(faces),
        covered,
    )
end

runs = [
    quotient_stats("theirs (hand)", pclf_theirs),
    quotient_stats("ours k=0", pclf_k0),
    quotient_stats("ours k=1", pclf_k1),
]

println()
@printf(
    "%-14s %8s %7s %9s %11s %10s %12s %9s\n",
    "certificate",
    "cells",
    "slices",
    "Σ facets",
    "max f/cell",
    "vol",
    "cells/vol",
    "build s",
)
for r in runs
    @printf(
        "%-14s %8d %7d %9d %11d %10.2f %12.2f %9.1f\n",
        r.name,
        r.ncells,
        r.nslices,
        r.Σfacets,
        r.maxfacets,
        r.covered,
        r.ncells / r.covered,
        r.t,
    )
end

# ---------------------------------------------------------
# Restricted to the common working set
# ---------------------------------------------------------
# Each certificate's slice family overshoots the working set differently — Γ_X is
# certificate-relative — so the whole-quotient counts above partly measure regions the problem
# never asked about. Every quotient covers X by construction, which makes X the canonical common
# domain: the counts below keep only the cells that meet it, and put every certificate's
# complexity over the same 139.24 units of volume.

function restricted_stats(name, quotient)
    Xh = UT._as_hpolytope(X)
    kept = [q for q in values(quotient.states) if !UT.is_disjoint(q.set, Xh)]
    faces = [PCQ.num_faces(q.set) for q in kept]
    return (; name, ncells = length(kept), Σfacets = sum(faces), maxfacets = maximum(faces))
end

volX = LazySets.volume(UT._as_hpolytope(X); backend = CDDLib.Library())

println()
println("restricted to the working set X (volume $(round(volX; digits = 2))):")
@printf(
    "%-14s %8s %9s %11s %12s\n",
    "certificate",
    "cells",
    "Σ facets",
    "max f/cell",
    "cells/vol",
)
for r in runs
    s = restricted_stats(r.name, r.quotient)
    @printf(
        "%-14s %8d %9d %11d %12.2f\n",
        s.name,
        s.ncells,
        s.Σfacets,
        s.maxfacets,
        s.ncells / volX,
    )
end
