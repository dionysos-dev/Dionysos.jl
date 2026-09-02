# How the quotient changes with the order of the De Bruijn graph.
#
# The same two-mode system and specification as `debruijn_graph.jl`, run at k = 0, 1, 2. For
# M = 2 modes a De Bruijn graph of order k has 2^k nodes, so this sweeps 1, 2 and 4 nodes: more
# nodes means a richer path-complete family, a tighter contraction rate, and a finer quotient.

include(joinpath(dirname(@__DIR__), "common.jl"))
using Spot

gr()

function run_debruijn(;
    k::Int = 1,
    dual::Bool = true,
    template_mode::Symbol = :rotation,
    θ::Float64 = π / 6,
    max_iter::Int = 1000,
    max_levels::Int = 100,
    max_slices::Int = 8,
    atol::Float64 = 1e-3,
    print_level::Int = 1,
    early_stop::Bool = false,
    N::Int = 50,
)
    (; f, problem, R1, R2) = two_mode_problem()

    graph = PCLF.generate_DeBruijn_edges(2, k; dual = dual)
    nodes = sort(collect(graph.verts))

    pclf = PCLF.compute_symmetric_2n_faces_polyhedral_pieces_pclf(
        f,
        graph,
        JuMP.optimizer_with_attributes(Clarabel.Optimizer, "max_iter" => max_iter);
        Gmats = rotation_templates(nodes; θ = θ, mode = template_mode),
        MLF = true,
        verbose = false,
    )

    println()
    println("------------------------------------------")
    println("k = $k, nodes = $nodes ($(length(nodes)) of them)")
    println("Computed JSR upper bound / contraction rate = ", pclf.JSRapprox)

    (; quotient, D) = build_quotient(
        problem,
        pclf;
        atol = atol,
        max_levels = max_levels,
        max_slices = max_slices,
        print_level = print_level,
    )

    display(plot_cells_by_node(quotient; nodes = nodes))

    φ = ltl"F(R1 & F(D))"
    x0 = SVector(2.3, 1.5)

    result = synthesize_cosafe_ltl(
        f,
        quotient,
        Dionysos.spot_stepper(φ),
        Dict(:D => D, :R1 => R1, :R2 => R2),
        Dict(:D => -1, :R1 => 1, :R2 => 2),
        x0;
        early_stop = early_stop,
        print_level = print_level,
        N = N,
    )

    display(plot_controllable_set(quotient, result.controllable_set, result.X, string(φ)))
    display(
        plot_partition_by_node(
            quotient,
            result.controllable_set,
            result.uncontrollable_set,
            result.X;
            nodes = nodes,
        ),
    )

    return (; quotient, result)
end

for k in 0:2
    run_debruijn(; k = k)
end
