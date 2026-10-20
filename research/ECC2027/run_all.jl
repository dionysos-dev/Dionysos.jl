# Runs all the experiments:  julia --project=. --threads=auto run_all.jl
# Environment variables: see common.jl.

include(joinpath(@__DIR__, "common.jl"))
# Start from an empty numbers.tex.
rm(joinpath(RESULTS, "numbers.tex"); force = true)
write_environment()

for script in ("E1_pendulum.jl", "E2_cartpole.jl")
    println("\n", "#"^84, "\n# ", script, "\n", "#"^84)
    include(joinpath(@__DIR__, "experiments", script))
end
println("\nresults are in ", RESULTS)
