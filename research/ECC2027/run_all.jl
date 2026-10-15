# Regenerate every number, table and figure of the paper.
#
#     julia --project=. --threads=auto run_all.jl
#
# Results go to `results/`; the paper quotes copies of them.
# Switches: SAMPLES=n for the median of n timed rounds (default 3), FIGURES=0 to skip the plots.

include(joinpath(@__DIR__, "common.jl"))
write_environment()

for script in ("E1_pendulum.jl", "E2_cartpole.jl")
    println("\n", "#"^84, "\n# ", script, "\n", "#"^84)
    include(joinpath(@__DIR__, "experiments", script))
end
println("\nresults are in ", RESULTS)
