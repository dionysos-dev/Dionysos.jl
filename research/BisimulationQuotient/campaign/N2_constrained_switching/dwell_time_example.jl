# Minimum dwell time: the best-motivated constrained-switching example.
#
# Two modes, each STRONGLY contracting on its own (rho = 0.3, they would settle in a few steps), whose
# fast alternation diverges (rho(A1*A2)^(1/2) = 1.95). So the system is unstable under arbitrary
# switching and no common Lyapunov function exists at any complexity -- the predecessor has no input.
# Impose a minimum dwell time and stability returns.
#
#     A1 = [0.3  2 ;  0   0.3]      A2 = [0.3  0 ; -2  0.3]
#
# Measured (conic orders 1/2/3):
#
#   unconstrained   1 node    2.04402  1.95393  1.95394     no certificate
#   dwell tau = 1   2 nodes   2.04403  1.95393  1.95392     constraint vacuous -- reproduces the above
#   dwell tau = 2   4 nodes   1.14971  1.11405  1.09234     still no certificate
#   dwell tau = 3   6 nodes   0.84874  0.83708  0.82515     *** certifies ***
#   dwell tau = 4   8 nodes   0.70558  0.70014  0.69458
#   dwell tau = 5  10 nodes   0.62271  0.61963  0.61651
#
# The threshold is VISIBLE: a dwell of 2 is not enough, a dwell of 3 is. "The minimum dwell time is
# three steps" is a number the experiment discovers rather than assumes, and the 6-node automaton is
# small enough to build a quotient on.
#
# The tau = 1 row is a free sanity check: the constraint is vacuous there and the rate reproduces the
# unconstrained one to six digits.
#
# Why dwell time rather than scaling a matrix: it is physically motivated (actuators cannot switch
# arbitrarily fast), it is the most classical language restriction there is (`PEDJ:16`, `ZHANG201342`,
# `Thesis:MPhilippe`), and it gives a knob with a threshold instead of a yes/no -- so the sweep shows
# how much conservatism must be accepted to buy a certificate, with the automaton growing as M*tau
# alongside.
#
# TRAP when retuning: `shear` must keep rho(A1*A2)^(1/2) > 1, or the system is not demonstrably
# unstable and the whole point collapses. A failing certificate on a stable system says only that the
# template is too weak, which is a much weaker claim.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

using LinearAlgebra

spectral_radius(M) = maximum(abs.(LinearAlgebra.eigvals(M)))

"""
Two shear modes, each contracting at `diag`, whose alternation diverges when `shear` is large enough.
"""
function shear_modes(; diag::Float64 = 0.3, shear::Float64 = 2.0)
    return [diag shear; 0.0 diag], [diag 0.0; -shear diag]
end

"""
The minimum-dwell-time automaton: state `(m, k)` is "in mode `m` for `k` steps, capped at `tau`".

Staying in the current mode is always allowed; switching out is allowed only from `k = tau`. The
generated language is exactly the words whose mode-runs have length at least `tau`, so `tau = 1`
imposes nothing and must reproduce the unconstrained answer.
"""
function dwell_graph(M::Int, tau::Int)
    edges = Tuple{Tuple{Int, Int}, Tuple{Int, Int}, Int}[]
    for m in 1:M, k in 1:tau
        push!(edges, ((m, k), (m, min(k + 1, tau)), m))
        if k == tau
            for m2 in 1:M
                m2 == m && continue
                push!(edges, ((m, k), (m2, 1), m2))
            end
        end
    end
    return PCLF.edgeList_to_LabDigraph(edges)
end

const OPTIMIZER = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 4000,
    "verbose" => false,
    "tol_feas" => 1e-6,
    "tol_gap_abs" => 1e-6,
    "tol_gap_rel" => 1e-6,
)

function rate(f, graph, order)
    partition = PCLF.conic_partitions_2d(order)
    value = try
        PCLF.compute_polyhedral_pieces_pclf(
            f,
            graph,
            OPTIMIZER,
            Dict(v => partition for v in graph.verts);
            MLF = true,
        ).JSRapprox
    catch
        Inf
    end
    return value
end

fmt(v) = isfinite(v) ? string(round(v; digits = 5)) : "none"

function report(; diag = 0.3, shear = 2.0, taus = 1:5, orders = 1:3)
    A1, A2 = shear_modes(; diag = diag, shear = shear)
    f = ST.with_switching(
        HybridSystems.discreteswitchedsystem([A1, A2]),
        HybridSystems.ControlledSwitching(),
    )

    println("A1 = ", A1, "   A2 = ", A2)
    println(
        "rho(A1) = ",
        round(spectral_radius(A1); digits = 4),
        "   rho(A2) = ",
        round(spectral_radius(A2); digits = 4),
        "   rho(A1*A2)^(1/2) = ",
        round(spectral_radius(A1 * A2)^0.5; digits = 4),
        spectral_radius(A1 * A2) > 1 ?
        "   <- alternation diverges: no common Lyapunov function exists" :
        "   <- NOT unstable; this example does not make the point",
    )
    println()

    header = [rpad("graph", 24), rpad("nodes", 7)]
    for o in orders
        push!(header, rpad("order $o", 10))
    end
    println(join(header))

    unconstrained = PCLF.edgeList_to_LabDigraph([(1, 1, 1), (1, 1, 2)])
    row = [rpad("unconstrained", 24), rpad(1, 7)]
    for o in orders
        push!(row, rpad(fmt(rate(f, unconstrained, o)), 10))
    end
    println(join(row))

    for tau in taus
        graph = dwell_graph(2, tau)
        row = [rpad("dwell tau = $tau", 24), rpad(length(graph.verts), 7)]
        for o in orders
            push!(row, rpad(fmt(rate(f, graph, o)), 10))
        end
        println(join(row))
    end
    return nothing
end

# Run the table only when this file is the entry point, so `plot_dwell_threshold.jl` can include it
# for `shear_modes`, `dwell_graph` and `rate` without paying for the sweep twice.
if abspath(PROGRAM_FILE) == @__FILE__
    report()
end
