# Runs all the experiments:  julia --project=. run_all.jl   [CACHE=0] [TIMING=0] [SAMPLES=n]
# Every script writes its numbers to results/numbers_<script>.tex, which the paper includes, and its
# figure to results/fig_<script>.png. Each runs in its own process: E2 runs a conference script whose
# constants would collide with those of common.jl.
scripts = (
    "E2_graph.jl",
    "E3_constrained_plant.jl",
    "E4_unstable_plant.jl",
    "cells_per_slice.jl",
    "E5_dwell_time.jl",
    "E6_coset_family.jl",
    "fig_coset_family.jl",
)
for script in scripts
    println("\n", "#"^84, "\n# ", script, "\n", "#"^84)
    run(
        `$(Base.julia_cmd()) --project=$(@__DIR__) $(joinpath(@__DIR__, "experiments", script))`,
    )
end
println("\nresults are in ", joinpath(@__DIR__, "results"))
