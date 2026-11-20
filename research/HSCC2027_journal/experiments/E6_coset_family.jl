# E6: the rotation family of Corollary 9.2 (paper, Sections 9 and 10.7).
#
#     julia --project=. experiments/E6_coset_family.jl   [TIMING=0] [SAMPLES=n]
#
# A_m = α U^(m-1), with U the rotation by π/(2M) in each of the n/2 planes. The coset certificate
# has M nodes and square pieces, rotated by s π/(2M) at node s, with rate α attained exactly; the
# induced common Lyapunov function is the product of regular 4M-gons. Both are written in closed form,
# without any LP. Without regions, Proposition 9.3 predicts M (N + 1) and N + 1 cells, the levels
# τ_0, …, τ_N that `N_FREE` counts; with regions, the two constructions are compared in cells,
# polytopes, facets and time, in dimensions 2 and 4.
include(joinpath(dirname(@__DIR__), "common.jl"))

const ALPHA = 0.8
const N_FREE = 5

rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

"The rotation family in dimension `n` with `M` modes and rate `α`, its coset certificate and its
common function."
function family(n, M; α = ALPHA)
    U = kron(Matrix(1.0LinearAlgebra.I, n ÷ 2, n ÷ 2), rotation(π / (2M)))
    As = [α * U^(m - 1) for m in 1:M]
    graph = PCLF.edgeList_to_LabDigraph([
        (s, mod(s + m - 1, M), m) for s in 0:(M - 1) for m in 1:M
    ])
    pieces = Dict{Any, PCLF.AbstractPiece}(
        s => PCLF.PolyhedralPiece(Matrix(transpose(U^s)), ones(n)) for s in 0:(M - 1)
    )
    pclf = PCLF.PCLF(graph, pieces, α)
    return (; As, α, pclf, common = PCLF.build_common_lyapunov(pclf; mode = :max))
end

box(low, high) = UT._as_hpolytope(LazySets.Hyperrectangle(; low, high))
"Two box regions away from the origin, padded with a slab in the extra coordinates."
function regions(n)
    pad(v, w) = [v; fill(w, n - 2)]
    return [
        box(pad([0.35, 0.10], -0.15), pad([0.65, 0.30], 0.15)),
        box(pad([-0.60, -0.55], -0.15), pad([-0.30, -0.25], 0.15)),
    ]
end

"The first level clearing every region, for both certificates, as in the other experiments."
function slices(fam, regs)
    clears(j) = all(
        c -> PCQ.all_nodes_clear_regions(c, fam.α^(j - 1), regs; tol = 1e-2),
        (fam.pclf, fam.common),
    )
    return something(findfirst(clears, 1:80), -1)
end

println("="^84, "\nE6 — the rotation family\n", "="^84)

# ── Proposition 9.3: no region, no split ────────────────────────────────────────────────────────────
closed_form = true
for M in (2, 3, 4)
    fam = family(2, M)
    f = switched(fam.As; controlled = true)
    problem = PR.BisimulationQuotientProblem(
        f,
        box(-ones(2), ones(2)),
        LazySets.HPolytope{Float64}[],
    )
    cert = (; N = N_FREE, ΓX = 1.0)
    ours = cell_stats(build("E6_free_M$(M)_pclf", problem, cert, fam.pclf).quotient)
    theirs = cell_stats(build("E6_free_M$(M)_common", problem, cert, fam.common).quotient)
    ok = ours.cells == M * N_FREE && theirs.cells == N_FREE
    global closed_form &= ok
    Printf.@printf(
        "no region, M = %d: %d cells on the graph (predicted %d), %d on the common function (predicted %d)  %s\n",
        M,
        ours.cells,
        M * N_FREE,
        theirs.cells,
        N_FREE,
        ok ? "ok" : "MISMATCH"
    )
end

# ── with regions, dimensions 2 and 4 ────────────────────────────────────────────────────────────────
rows = []
# In dimension 4, every set operation is far more expensive, and the rate 0.8 needs N = 6; with 0.6
# the regions are cleared with fewer slices, and the comparison stays within reach.
const ALPHA_4D = 0.6
for (n, M, α) in ((2, 2, ALPHA), (2, 3, ALPHA), (2, 4, ALPHA), (4, 2, ALPHA_4D))
    fam = family(n, M; α)
    f = switched(fam.As; controlled = true)
    regs = regions(n)
    problem = PR.BisimulationQuotientProblem(f, box(-ones(n), ones(n)), regs)
    cert = (; N = slices(fam, regs), ΓX = 1.0)
    tag = α == ALPHA ? "" : "_a$(round(Int, 100α))"
    built_common = build("E6_n$(n)_M$(M)$(tag)_common", problem, cert, fam.common)
    built_pclf = build("E6_n$(n)_M$(M)$(tag)_pclf", problem, cert, fam.pclf)
    push!(
        rows,
        (;
            n,
            M,
            α,
            N = cert.N,
            problem,
            cert,
            fam,
            common = cell_stats(built_common.quotient),
            pclf = cell_stats(built_pclf.quotient),
        ),
    )
    Printf.@printf(
        "n = %d, M = %d, α = %.1f, %d slices: common %d cells / %d polytopes / %.1f mean facets / %d max; graph %d / %d / %.1f / %d\n",
        n,
        M,
        α,
        cert.N,
        rows[end].common.cells,
        rows[end].common.polytopes,
        rows[end].common.mean_facets,
        rows[end].common.max_facets,
        rows[end].pclf.cells,
        rows[end].pclf.polytopes,
        rows[end].pclf.mean_facets,
        rows[end].pclf.max_facets
    )
end
# Timed once every quotient exists.
times = [
    (
        build_seconds(r.problem, r.cert, r.fam.common),
        build_seconds(r.problem, r.cert, r.fam.pclf),
    ) for r in rows
]

open(joinpath(RESULTS, "table_E6_coset_family.tex"), "w") do io
    println(io, "% Written by experiments/E6_coset_family.jl. Do not edit.")
    println(io, "\\begin{tabular}{ccc c rr rr rr rr}")
    println(io, "\\toprule")
    println(
        io,
        "& & & & \\multicolumn{2}{c}{cells} & \\multicolumn{2}{c}{polytopes} & \\multicolumn{2}{c}{max facets per cell} & \\multicolumn{2}{c}{time (s)} \\\\",
    )
    println(
        io,
        "\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}",
    )
    println(
        io,
        "\$n\$ & \$M\$ & \$\\alpha\$ & \$N\$ & iCLF & PCLF & iCLF & PCLF & iCLF & PCLF & iCLF & PCLF \\\\",
    )
    println(io, "\\midrule")
    for (r, (tc, tp)) in zip(rows, times)
        Printf.@printf(
            io,
            "%d & %d & %.1f & %d & %d & %d & %d & %d & %d & %d & %s & %s \\\\\n",
            r.n,
            r.M,
            r.α,
            r.N - 1,
            r.common.cells,
            r.pclf.cells,
            r.common.polytopes,
            r.pclf.polytopes,
            r.common.max_facets,
            r.pclf.max_facets,
            fmt_seconds(tc),
            fmt_seconds(tp)
        )
    end
    println(io, "\\bottomrule")
    return println(io, "\\end{tabular}")
end
write_numbers(
    "E6_coset_family",
    Pair{String, Any}[
        "ESixClosedFormHolds" => closed_form ? "true" : "false",
        "ESixAlpha" => fmt(ALPHA),
        "ESixFreeN" => fmt(N_FREE - 1),
    ],
)

# The illustration of the planar case M = 2 is drawn by experiments/fig_coset_family.jl.
println("table and numbers written to results/")
