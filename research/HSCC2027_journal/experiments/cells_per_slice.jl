# The number of cells in each slice of the quotients of E3 and E4 (paper, Section 10.5).
#
#     julia --project=. experiments/cells_per_slice.jl
#
# Reads the cached quotients, so E3 and E4 must have run first. The bound of Theorem 8.2 on the cells
# of a slice grows by a factor tending to λⁿ per slice, with λ the growth rate of the admissible
# language: 2 for arbitrary switching and (1+√5)/2 under "no 11". The counts lie far below the bound,
# so their own factor may exceed λⁿ in the inner slices; it falls below it in the outer ones.
include(joinpath(dirname(@__DIR__), "common.jl"))

quotients = [
    ("arbitrary switching", "Free", "E3_free_6.jld2", :black),
    ("constrained plant", "Cons", "E3_constrained_6.jld2", WON),
    ("restricted unstable plant", "Unstable", "E4_restricted_order3_atol4_9.jld2", MISSED),
]

figure = Plots.plot(;
    yscale = :log10,
    xlabel = "slice index",
    ylabel = "cells in the slice",
    legend = :topleft,
    framestyle = :box,
    size = (560, 400),
    dpi = 200,
)
numbers = Pair{String, Any}[]
for (label, tag, file, colour) in quotients
    path = joinpath(CACHE_DIR, file)
    isfile(path) || error("missing $file: run E3 and E4 first")
    quotient = JLD2.jldopen(path, "r") do io
        return io["quotient"]
    end
    # The code numbers the slices from 1, the paper from 0.
    slices = [state.slice for state in values(quotient.states)]
    counts = [count(==(i), slices) for i in 1:maximum(slices)]
    # The factor between consecutive slices: its peak, and its value between the last two slices.
    growth = counts[2:end] ./ counts[1:(end - 1)]
    Printf.@printf(
        "%-28s cells per slice %s, growth factors %s\n",
        label,
        counts,
        join([Printf.@sprintf("%.2f", g) for g in growth], " ")
    )
    Plots.plot!(
        figure,
        0:(length(counts) - 1),
        counts;
        label,
        color = colour,
        marker = :circle,
        linewidth = 1.6,
    )
    push!(numbers, "SliceGrowthPeak$tag" => fmt(maximum(growth)))
    push!(numbers, "SliceGrowth$tag" => fmt(growth[end]))
end
Plots.savefig(figure, joinpath(RESULTS, "fig_cells_per_slice.png"))
write_numbers("cells_per_slice", numbers)
println("figure written to results/fig_cells_per_slice.png")
