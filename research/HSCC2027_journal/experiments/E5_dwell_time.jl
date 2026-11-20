# E5: two modes that contract on their own and diverge when alternated, under a minimum dwell time
# (paper, Section 10.6).
#
#     julia --project=. experiments/E5_dwell_time.jl
#
# The dwell time belongs to the plant, and its constraint graph represents itself. For each dwell
# time τ, the certified rate of Proposition 3.1 on that graph, for three templates, is compared with
# a lower bound on the constrained joint spectral radius from the cycles 1^a 2^b of the graph. The
# smallest τ that certifies is the threshold; the lower bound exceeding one below it shows that the
# threshold belongs to the system and not to the templates.
include(joinpath(dirname(@__DIR__), "common.jl"))

const A_DWELL = [[0.3 2.0; 0.0 0.3], [0.3 0.0; -2.0 0.3]]
const TAUS = 1:6
const ORDERS = 1:3

f = switched(A_DWELL; controlled = false)

"Minimum dwell time `τ` over two modes: node `(m, k)` records that mode `m` has been applied `k` times
in a row, capped at `τ`, and another mode may follow only at `k = τ`."
function dwell_graph(τ)
    edges = Tuple{Tuple{Int, Int}, Tuple{Int, Int}, Int}[]
    for m in 1:2, k in 1:τ
        push!(edges, ((m, k), (m, min(k + 1, τ)), m))
        k == τ && push!(edges, ((m, τ), (3 - m, 1), 3 - m))
    end
    return PCLF.edgeList_to_LabDigraph(edges)
end

spectral(M) = maximum(abs, LinearAlgebra.eigvals(M))
"Lower bound on the constrained joint spectral radius: the cycles 1^a 2^b with τ ≤ a, b ≤ τ + 8."
lower_bound(τ) = maximum(
    spectral(A_DWELL[2]^b * A_DWELL[1]^a)^(1 / (a + b)) for a in τ:(τ + 8) for
    b in τ:(τ + 8)
)

"The certified rate on `graph` with the uniform conic template of `order`, `Inf` when the LP fails."
function certified_rate(graph, order)
    cones = PCLFProducer.conic_partitions_2d(order)
    pclf = try
        PCLFProducer.compute_polyhedral_pieces_pclf(
            f,
            graph,
            HiGHS.Optimizer,
            Dict(nd => cones for nd in graph.verts);
            MLF = true,
        )
    catch
        return Inf
    end
    return pclf.JSRapprox
end

println("="^84, "\nE5 — minimum dwell time\n", "="^84)
alternation = spectral(A_DWELL[1] * A_DWELL[2])^(1 / 2)
Printf.@printf(
    "ρ(A1) = ρ(A2) = %.2f, ρ(A1 A2)^(1/2) = %.4f\n",
    spectral(A_DWELL[1]),
    alternation
)
@assert PCLF.respects_constraint(
    dwell_graph(3),
    Dict((m, 3) => [(m, 3)] for m in 1:2),
    dwell_graph(3),
)

rates = [certified_rate(dwell_graph(τ), order) for τ in TAUS, order in ORDERS]
bounds = [lower_bound(τ) for τ in TAUS]
for (i, τ) in enumerate(TAUS)
    Printf.@printf(
        "τ = %d  (%2d nodes)  lower bound %.4f   rates %s\n",
        τ,
        2τ,
        bounds[i],
        join([Printf.@sprintf("%.4f", r) for r in rates[i, :]], "  ")
    )
end
best = vec(minimum(rates; dims = 2))
threshold = TAUS[findfirst(<(1), best)]
Printf.@printf(
    "threshold: τ* = %d, best rate %.4f; lower bound at τ* - 1: %.4f\n",
    threshold,
    best[threshold],
    threshold > 1 ? bounds[threshold - 1] : NaN
)

figure = Plots.plot(;
    xlabel = "minimum dwell time τ",
    ylabel = "rate",
    legend = :topright,
    framestyle = :box,
    size = (560, 400),
    dpi = 200,
)
Plots.hline!(figure, [1.0]; color = :gray, linestyle = :dot, label = false)
Plots.plot!(
    figure,
    TAUS,
    bounds;
    color = :black,
    linestyle = :dash,
    linewidth = 1.4,
    label = "lower bound (cycles)",
)
for (j, order) in enumerate(ORDERS)
    Plots.plot!(
        figure,
        TAUS,
        rates[:, j];
        marker = :circle,
        linewidth = 1.4,
        label = "certified, template of order $order",
    )
end
Plots.savefig(figure, joinpath(RESULTS, "fig_E5_dwell_time.png"))

write_numbers(
    "E5_dwell_time",
    Pair{String, Any}[
        "EFiveAlternation" => fmt(alternation; digits = 3),
        "EFiveThreshold" => fmt(threshold),
        "EFiveRateStar" => fmt(best[threshold]; digits = 3),
        "EFiveBoundBefore" => fmt(bounds[threshold - 1]; digits = 3),
        "EFiveBoundStar" => fmt(bounds[threshold]; digits = 3),
    ],
)
println(
    "numbers written to results/numbers_E5_dwell_time.tex, figure to results/fig_E5_dwell_time.png",
)
