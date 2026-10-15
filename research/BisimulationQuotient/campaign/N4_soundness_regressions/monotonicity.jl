# Is the certified rate non-increasing over nested templates?
#
# Conic partition orders 1 ⊂ 2 ⊂ 3 give nested templates, so the rate a certificate achieves must
# be non-increasing in the order. A violation means the conic program is returning something that
# is not the optimum, and that is exactly the artefact that forced a retraction once before.
#
# Measured result: monotone in all six cases at the loosened tolerances, and non-monotone in five
# of six at Clarabel's defaults — with the failure confined to order 3, where the default 1e-8
# feasibility tolerance stops converging and returns a rate worse than order 1, or `Inf`. Orders 1
# and 2 agree to the fifth decimal under both settings, so no published number is affected.
#
# The conclusion is a solver requirement, not a bug: any experiment at conic order ≥ 3 must set
# `tol_feas`/`tol_gap_*` to 1e-6. This file is the regression test for it.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

const LOOSE = JuMP.optimizer_with_attributes(
    Clarabel.Optimizer,
    "max_iter" => 4000,
    "tol_feas" => 1e-6,
    "tol_gap_abs" => 1e-6,
    "tol_gap_rel" => 1e-6,
)

const DEFAULTS = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 4000)

"""
Rates certified over conic partition orders 1, 2, 3 by the De Bruijn graph of order `k`.

`Inf` marks a certificate the solver failed to find, which is how the solver reports
infeasibility rather than returning a bound.
"""
function rates_over_orders(f, k, optimizer)
    graph = PCLF.generate_DeBruijn_edges(2, k)
    out = Float64[]
    for order in 1:3
        partition = PCLF.conic_partitions_2d(order)
        rate = try
            PCLF.compute_polyhedral_pieces_pclf(
                f,
                graph,
                optimizer,
                Dict(v => partition for v in graph.verts);
                MLF = true,
            ).JSRapprox
        catch
            Inf
        end
        push!(out, rate)
    end
    return out
end

is_monotone(r) = all(r[i] >= r[i + 1] - 1e-9 for i in 1:(length(r) - 1))

fmt_rate(v::Float64) = isfinite(v) ? string(round(v; digits = 6)) : "Inf"

function report(systems, tolerances)
    for (tolname, optimizer) in tolerances
        println("\n===== tolerance: ", tolname, " =====")
        println(
            rpad("system", 18),
            rpad("nodes", 7),
            rpad("order1", 11),
            rpad("order2", 11),
            rpad("order3", 11),
            "monotone?",
        )
        for (name, f) in systems, k in 0:1
            r = rates_over_orders(f, k, optimizer)
            println(
                rpad(name, 18),
                rpad(2^k, 7),
                rpad(fmt_rate(r[1]), 11),
                rpad(fmt_rate(r[2]), 11),
                rpad(fmt_rate(r[3]), 11),
                is_monotone(r) ? "yes" : "*** NO ***",
            )
        end
    end
    return nothing
end

systems = [
    ("gol_lazar_belta", gol_lazar_belta_problem().f),
    ("observer_graph", observer_graph_problem(; p = 1.7).f),
    ("two_mode", two_mode_problem().f),
]

report(systems, [("loosened 1e-6", LOOSE), ("Clarabel defaults", DEFAULTS)])
