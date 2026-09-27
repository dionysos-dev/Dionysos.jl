# Run the three experiments end to end. This is the entry point for artifact evaluation.
#
#     julia --project=. run_all.jl
#
# Expected total runtime: a few minutes on a laptop. Set SAMPLES=3 for medians over three timed
# rounds instead of one, and FIGURES=0 to skip every figure.

for script in (
    "experiment1_gol_lazar_belta.jl",
    "experiment2_graph_orientation.jl",
    "experiment3_diverse_pieces.jl",
)
    println("\n", "#"^84)
    println("# ", script)
    println("#"^84)
    include(joinpath(@__DIR__, script))
end
println("\nfigures are in ", joinpath(@__DIR__, "figures"))
