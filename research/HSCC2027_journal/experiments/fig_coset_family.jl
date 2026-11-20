# The rotation family of Corollary 9.2 for n = 2 and M = 2 (paper, Section 10.7): the square piece of
# node 0, the diamond piece of node 1, and the octagon of the common Lyapunov function they induce.
#
#     julia --project=. experiments/fig_coset_family.jl
include(joinpath(dirname(@__DIR__), "common.jl"))

U = [cos(π / 4) -sin(π / 4); sin(π / 4) cos(π / 4)]
square = PCLF.get_sublevel_set(
    PCLF.PolyhedralPiece(Matrix(1.0LinearAlgebra.I, 2, 2), ones(2)),
    1.0,
)
diamond = PCLF.get_sublevel_set(PCLF.PolyhedralPiece(Matrix(transpose(U)), ones(2)), 1.0)
octagon = LazySets.intersection(square, diamond)

# No legend: the caption names the three shapes, which their outlines tell apart.
figure = Plots.plot(;
    aspect_ratio = :equal,
    framestyle = :box,
    legend = false,
    size = (440, 440),
    dpi = 200,
    xlims = (-1.55, 1.55),
    ylims = (-1.55, 1.55),
)
Plots.plot!(
    figure,
    outline(octagon)...;
    seriestype = :shape,
    fillcolor = :gray85,
    linecolor = :black,
    linewidth = 1.8,
)
Plots.plot!(figure, outline(square)...; color = WON, linestyle = :dash, linewidth = 1.6)
Plots.plot!(figure, outline(diamond)...; color = MISSED, linestyle = :dash, linewidth = 1.6)
Plots.savefig(figure, joinpath(RESULTS, "fig_E6_coset_family.png"))
println("figure written to results/fig_E6_coset_family.png")
