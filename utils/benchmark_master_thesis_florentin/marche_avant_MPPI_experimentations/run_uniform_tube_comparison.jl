# Thin entry point kept next to the original marche-avant benchmark.
# The implementation lives in ../corridor_counter_example so the new scripts
# and generated artifacts stay grouped in the corridor-counter-example folder.

include(
    joinpath(@__DIR__, "..", "corridor_counter_example", "run_uniform_tube_comparison.jl"),
)

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
