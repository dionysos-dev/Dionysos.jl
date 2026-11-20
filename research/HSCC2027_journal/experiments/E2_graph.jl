# E2: the primal and the dual graph of order one, with almost identical pieces (paper, Section 10.2).
#
#     julia --project=. experiments/E2_graph.jl
#
# Runs the conference experiment unchanged, research/HSCC2027/experiment2_primal_and_dual.jl, with
# SAMPLES=3, which its README requires to separate the build times of its shortest builds, and
# writes its counts and times for the paper.
ENV["SAMPLES"] = get(ENV, "SAMPLES", "3")
ENV["FIGURES"] = "0"
include(joinpath(dirname(dirname(@__DIR__)), "HSCC2027", "experiment2_primal_and_dual.jl"))

"A build time with three significant digits, as `fmt_seconds` of the capsule's common.jl."
seconds(t) = isnan(t) ? "--" : @sprintf("%.*f", max(0, 2 - floor(Int, log10(t))), t)

open(joinpath(dirname(@__DIR__), "results", "numbers_E2_graph.tex"), "w") do io
    println(io, "% Written by experiments/E2_graph.jl. Do not edit.")
    for (r, graph) in ((primal, "Primal"), (dual, "Dual")),
        (a, arm) in ((r.a1, "Common"), (r.a2, "Graph"))

        name = "ETwo$(graph)$(arm)"
        println(io, "\\newcommand{\\$(name)Cells}{$(a.s.cells)}")
        println(io, "\\newcommand{\\$(name)Polytopes}{$(sum(a.s.parts))}")
        println(io, "\\newcommand{\\$(name)MeanFacets}{$(@sprintf("%.1f", a.s.mean_f))}")
        println(io, "\\newcommand{\\$(name)MaxFacets}{$(a.s.max_f)}")
        println(io, "\\newcommand{\\$(name)Time}{$(seconds(a.t))}")
    end
    # The cells of the induced function against those of one node of the dual graph, which has two.
    factor = dual.a1.s.cells / (dual.a2.s.cells / 2)
    return println(io, "\\newcommand{\\ETwoDualNodeFactor}{$(@sprintf("%.1f", factor))}")
end
println("numbers written to results/numbers_E2_graph.tex")
