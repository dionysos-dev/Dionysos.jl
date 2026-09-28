# What does route 1 actually pay to induce the common Lyapunov function?
#
# The cost splits in two, and both halves are measured here.
#
# The determinization is free, and the 2^|S| bound is nowhere near attained:
#
#   E2, 4-node PCLF   build_common_lyapunov  96 µs    observer states  3 of 16
#   E8, 2-node PCLF   build_common_lyapunov  86 µs    observer states  3 of 4
#
# E2's observer collapses to {1,2,3,4} -> {1,3,4} -> {2,4}: *fewer* states than the graph has nodes.
# So the exponential determinization the cost argument leans on is a worst case these certificates
# do not approach, and it should be presented as such.
#
# The geometry is where route 1 genuinely pays. `get_sublevel_set(::ObserverCLFPiece, γ)` builds
# ⋃_Q ⋂_{s∈Q} P_Γ^(s) by repeated `set_difference_decompose` with a `clean_poly` LP pass on each
# intersection, where route 2 emits an H-polytope directly with no LP at all: 150x to 1000x per
# sublevel set, a ratio far exceeding both the facet ratio and the piece ratio, which is the
# multiplicative-piece mechanism showing up exactly where predicted.
#
# And yet it does not matter at these sizes: sublevel sets are built once per level, so at 20-50
# levels the extra cost is 0.3-1.2 s against route 1's measured 41.5 s and 168.7 s -- under 1 %.
# Including the induced-common construction in route 1's column is correct and changes nothing.
#
# `th:induced_common` case 1 would avoid most of this on a COMPLETE graph, where V* = min_s V_s. On
# E8 the graph is complete and the pieces are nested, so that is a single 8-row polytope against the
# 21 pieces the general path builds -- route 1 is being measured while handicapped. Adding the
# shortcut is what this experiment is for.

include(joinpath(dirname(dirname(@__DIR__)), "common.jl"))

n_facets(P) = length(LazySets.constraints_list(P))

parts_of(S::LazySets.UnionSetArray) = S.array
parts_of(P::LazySets.HPolytope) = [P]

"""
Determinization cost and observer size, then the per-level geometry cost of both routes.
"""
function report(name, pclf)
    println("\n===== ", name, " =====")
    println("nodes = ", length(pclf.pieces), "   rate = ", pclf.JSRapprox)

    PCLF.build_common_lyapunov(pclf)   # warm-up, so the timing is not compilation
    t_common = @elapsed clf = PCLF.build_common_lyapunov(pclf)
    observer = clf.pieces[:clf].observer_states
    println(
        "build_common_lyapunov: ",
        round(t_common * 1e6; digits = 1),
        " µs    observer states = ",
        length(observer),
        " of max 2^",
        length(pclf.pieces),
        " = ",
        2^length(pclf.pieces),
    )

    for γ in (1.0, 0.5)
        t_nodes, n_nodes, f_nodes = 0.0, 0, 0
        for (_, piece) in pclf.pieces
            PCLF.get_sublevel_set(piece, γ)
            t_nodes += @elapsed P = PCLF.get_sublevel_set(piece, γ)
            parts = parts_of(P)
            n_nodes += length(parts)
            f_nodes += sum(n_facets, parts)
        end

        PCLF.get_sublevel_set(clf.pieces[:clf], γ)
        t_common_geom = @elapsed C = PCLF.get_sublevel_set(clf.pieces[:clf], γ)
        parts = parts_of(C)

        println(
            "  γ=",
            γ,
            " | route 2 (all nodes): ",
            round(t_nodes * 1e3; digits = 2),
            " ms, ",
            n_nodes,
            " pieces, ",
            f_nodes,
            " facets",
            " || route 1 (induced): ",
            round(t_common_geom * 1e3; digits = 2),
            " ms, ",
            length(parts),
            " pieces, ",
            sum(n_facets, parts),
            " facets   => time x",
            round(t_common_geom / max(t_nodes, 1e-9); digits = 1),
        )
    end
    return nothing
end

(; f) = observer_graph_problem(; p = 1.7)
report("E2  observer-graph problem, 4-node PCLF", observer_graph_pclf(f))

(; f) = gol_lazar_belta_problem()
report(
    "E8  Gol-Lazar-Belta problem, 2-node PCLF",
    PCLF.compute_polyhedral_pieces_pclf(
        f,
        PCLF.generate_DeBruijn_edges(2, 1),
        JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => 1000),
        Dict((1,) => PCLF.conic_partitions_2d(2), (2,) => PCLF.conic_partitions_2d(2));
        MLF = true,
    ),
)
