# Run the three experiments end to end. This is the entry point for artifact evaluation.
#
#     julia --project=. run_all.jl
#
# About twenty minutes from a cold `cache/`, most of it experiment 1. The switches are the ones the
# experiments take: SAMPLES=3 for medians over three timed rounds instead of one, FIGURES=0 to skip
# every figure, TIMING=0 to keep the counts and drop the measurement, CACHE=0 to ignore `cache/`.
#
# The figures this writes are GR's, in Computer Modern. The paper's are the same plots typeset by
# LaTeX: run `paper_figures/paper_figures_exp{1,2,3}.jl` afterwards, which needs `pdflatex` on PATH.

for script in (
    "experiment1_gol_lazar_belta.jl",
    "experiment2_primal_and_dual.jl",
    "experiment3_diverse_pieces.jl",
)
    println("\n", "#"^84)
    println("# ", script)
    println("#"^84)
    include(joinpath(@__DIR__, script))
end
println("\nfigures are in ", joinpath(@__DIR__, "figures"))
