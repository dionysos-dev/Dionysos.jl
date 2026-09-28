# The separation does not need switching: one mode, and the graph node is a clock.
#
#     x⁺ = A x,   A = ρ R(2π/q)          -- a single linear system, no switching
#     graph  the q-cycle  s -> s+1  over the one-letter alphabet, path-complete for 1*
#     PCLF   V_s(x) = ‖R^{-s} x‖_∞        -- four facets per node, for every q
#
# Along the edge (s, s+1),
#     V_{s+1}(Ax) = ‖R^{-(s+1)} ρ R x‖_∞ = ρ ‖R^{-s} x‖_∞ = ρ V_s(x),
# so the cycle certifies at rate exactly ρ. A one-step polytopic Lyapunov function needs
# R B ⊆ (γ/ρ) B; as γ → ρ this forces B invariant under the order-q rotation group, and a polytopic
# B then needs at least q facets. Restricted to four, the rate is bounded away from ρ.
#
# That is `memory_vs_geometry.jl`'s orbit argument with the second mode deleted -- it never used it.
#
# Measured result: the single node is above 1 for every (q, ρ) tried, so it certifies *nothing*,
# while the cycle attains ρ exactly. The ratio single/ρ reproduces the analytic gauge constants:
# 1.2599 at q = 5, saturating at 1.4141 ≈ √2 by q = 8 (a square rotated by 45°).
#
# With M = 1 the node is a clock: the memory is the phase k mod q, so the method becomes a periodic
# Lyapunov function abstraction. The certificate decreases every q steps while the abstraction still
# respects per-step transitions and per-step observations -- which is why the q steps cannot simply
# be lumped into one.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using Random
using LinearAlgebra

rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

# One mode only. Not a switched system.
function one_mode_system(ρ, q)
    return HybridSystems.discreteswitchedsystem([ρ .* rotation(2π / q)])
end

"""
The `q`-cycle over the single letter, with the analytic pieces `V_s(x) = ‖R^{-s}x‖_∞`.
"""
function cycle_pclf(ρ, q)
    θ = 2π / q
    graph = PCLF.edgeList_to_LabDigraph([(s, mod(s + 1, q), 1) for s in 0:(q - 1)])
    pieces = Dict{Any, PCLF.AbstractPiece}(
        s => PCLF.PolyhedralPiece(Matrix(rotation(-s * θ)), ones(2)) for s in 0:(q - 1)
    )
    return PCLF.PCLF(graph, pieces, ρ)
end

"""
Largest observed `V_{s+1}(Ax) / V_s(x)` over the cycle, sampled at random points.

An independent check on the analytic certificate involving no solver: it must return `ρ`.
"""
function verify_cycle_rate(ρ, q; nsample = 4000, seed = 3)
    θ = 2π / q
    A = ρ .* rotation(2π / q)
    V(s, x) = maximum(abs.(rotation(-s * θ) * x))
    rng = Random.MersenneTwister(seed)
    worst = 0.0
    for _ in 1:nsample
        x = randn(rng, 2)
        for s in 0:(q - 1)
            vs = V(s, x)
            vs < 1e-12 && continue
            worst = max(worst, V(mod(s + 1, q), A * x) / vs)
        end
    end
    return worst
end

"""
Best rate a single node certifies with a four-facet template, over `draws` orientations.

The identity is always included: a purely random search can miss the axis-aligned template and
would then understate what one node achieves, which is how a separation gets manufactured.
"""
function single_node_rate(ρ, q; draws = 40, seed = 11)
    f = one_mode_system(ρ, q)
    graph = PCLF.edgeList_to_LabDigraph([(1, 1, 1)])
    optimizer = JuMP.optimizer_with_attributes(
        Clarabel.Optimizer,
        "max_iter" => 3000,
        "verbose" => false,
        "tol_feas" => 1e-6,
        "tol_gap_abs" => 1e-6,
        "tol_gap_rel" => 1e-6,
    )
    rng = Random.MersenneTwister(seed)
    best = Inf
    for i in 0:draws
        M = i == 0 ? Matrix{Float64}(I, 2, 2) : Matrix(qr(randn(rng, 2, 2)).Q)
        rate = try
            PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
                f,
                graph,
                optimizer;
                Gmats = Dict(v => M for v in graph.verts),
                MLF = true,
                verbose = false,
            ).JSRapprox
        catch
            Inf
        end
        best = min(best, rate)
    end
    return best
end

function report(qs, ρs)
    println("ONE mode: x⁺ = ρ R(2π/q) x.  Four facets per node throughout.")
    println(
        rpad("q", 4),
        rpad("ρ", 7),
        rpad("1 node", 12),
        rpad("q-cycle", 10),
        rpad("verified", 11),
        "separation?",
    )
    for q in qs, ρ in ρs
        one = single_node_rate(ρ, q)
        checked = verify_cycle_rate(ρ, q)
        separates = (!isfinite(one) || one >= 1.0) && ρ < 1.0
        println(
            rpad(q, 4),
            rpad(ρ, 7),
            rpad(isfinite(one) ? string(round(one; digits = 5)) : "none", 12),
            rpad(round(ρ; digits = 5), 10),
            rpad(round(checked; digits = 5), 11),
            separates ? "*** YES ***" : "no",
        )
    end
    return nothing
end

# Guarded so `2_quotient_on_clock.jl` can include this file for `one_mode_system`, `cycle_pclf` and
# `single_node_rate` without paying for the 41-draw sweep at every q and rho.
if abspath(PROGRAM_FILE) == @__FILE__
    report([5, 6, 7, 8], [0.85, 0.90, 0.95])
end
