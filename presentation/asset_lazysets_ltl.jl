# Presentation asset: generic LazySets all the way through a co-safe LTL task.
# The hazard is a lazy union (zonotope ∪ box), the waypoint and goal are balls;
# the formula asks to visit the waypoint, then the goal, never touching the
# hazard. Output: presentation/gallery/lazysets_ltl.png
#
# Run: julia --project=test presentation/asset_lazysets_ltl.jl

ENV["GKSwstype"] = "100"

using StaticArrays, JuMP, Plots
import LazySets
using Symbolics, MathOptSymbolicAD
using Spot
using Dionysos
const MP = Dionysos.Mapping
const UT = Dionysos.Utils

model = Model(Dionysos.Optimizer)
@variable(model, -5.0 <= x1 <= 5.0)
@variable(model, -5.0 <= x2 <= 5.0)
@variable(model, -1.0 <= u1 <= 1.0)
@variable(model, -1.0 <= u2 <= 1.0)
@constraint(model, ∂(x1) == u1)
@constraint(model, ∂(x2) == u2)
@constraint(model, start(x1) in MOI.Interval(-3.6, -3.4))
@constraint(model, start(x2) in MOI.Interval(-3.6, -3.4))

hazard_set = LazySets.UnionSetArray([
    LazySets.Zonotope([-0.5, 0.0], [0.6 1.1; 1.3 -0.5]),
    LazySets.Hyperrectangle(; low = [1.0, -3.2], high = [2.2, -0.8]),
])
# a tilted ellipsoid (shape matrix Q)
waypoint_set = LazySets.Ellipsoid([-3.0, 3.2], [0.9 0.5; 0.5 0.6])
# a non-convex goal: a box with a ball bitten out of it (lazy set difference)
goal_box = LazySets.Hyperrectangle(; low = [2.4, 2.2], high = [4.2, 4.2])
goal_bite = LazySets.Ball2([2.4, 3.2], 0.9)
goal_set = UT.set_minus(goal_box, goal_bite)

@constraint(model, wp, [x1, x2] in Label(waypoint_set))
@constraint(model, gl, [x1, x2] in Label(goal_set))
@constraint(model, hz, [x1, x2] in Label(hazard_set; semantics = MP.OUTER))
@specification(model, ltl"F(wp & F(gl)) & G(!hz)")

set_attribute(model, "jacobian_bound", u -> SMatrix{2, 2}(0.0, 0.0, 0.0, 0.0))
set_attribute(model, "time_step", 0.3)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.2, 0.2)))
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))
set_attribute(model, "print_level", 0)
optimize!(model)
println("termination: ", termination_status(model))

trajectory = Dionysos.simulate(model, SVector(-3.5, -3.5); nsteps = 200)

fig = plot(; aspect_ratio = :equal, legend = false, dpi = 200)
plot!(
    fig,
    LazySets.Hyperrectangle(; low = [-5.0, -5.0], high = [5.0, 5.0]);
    color = :lightgray,
    linecolor = :gray,
)
plot!(fig, hazard_set; color = :black, opacity = 0.8)
plot!(fig, waypoint_set; color = :orange, opacity = 0.6)
# draw the set difference as box-then-bite (the lazy Intersection has no recipe)
plot!(fig, goal_box; color = :green, opacity = 0.5)
plot!(fig, goal_bite; color = :lightgray, opacity = 1.0, linecolor = :gray)
plot!(fig, trajectory; ms = 1.2, color = :blue)
annotate!(fig, -3.0, 4.45, text("wp", 11))
annotate!(fig, 3.6, 4.55, text("gl", 11))
annotate!(fig, -0.5, -2.2, text("hz", 11))
out = joinpath(@__DIR__, "gallery", "lazysets_ltl.png")
savefig(fig, out)
println("saved ", out)
