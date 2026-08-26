# The quotient over a state space that is not axis-aligned.
#
# Same system and observer graph as `pclf_vs_clf.jl`, but the state-space polytope is rotated by
# 10°, so the sublevel sets no longer line up with it. The figure shows the resulting cells over
# the problem geometry.

include(joinpath(@__DIR__, "common.jl"))

gr()

(; f, problem) = observer_graph_problem(; p = 2.5, rotation = deg2rad(10.0))

pclf_poly = observer_graph_pclf(f)
println("Computed JSR upper bound / contraction rate = ", pclf_poly.JSRapprox)

(; quotient) =
    build_quotient(problem, pclf_poly; atol = 1e-4, level_tol = 1e-2, max_slices = 30)

fig = plot(; aspect_ratio = :equal)
plot!(fig, quotient; what = :states, node = 1, show_contours = false)
plot!(fig, problem; opacity = 0.2)
display(fig)
