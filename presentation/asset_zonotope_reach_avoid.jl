# Presentation asset: reach-avoid straight through the JuMP front-end with generic
# LazySets — a zonotope obstacle and a ball target, no conversion, no special casing.
# Output: presentation/zonotope_reach_avoid.png
#
# Run: julia --project=test presentation/asset_zonotope_reach_avoid.jl

ENV["GKSwstype"] = "100"

using StaticArrays, JuMP, Plots
import LazySets
using Symbolics, MathOptSymbolicAD
using Dionysos
const MP = Dionysos.Mapping

model = Model(Dionysos.Optimizer)
@variable(model, -5.0 <= x1 <= 5.0)
@variable(model, -5.0 <= x2 <= 5.0)
@variable(model, -1.0 <= u1 <= 1.0)
@variable(model, -1.0 <= u2 <= 1.0)
@constraint(model, ∂(x1) == u1)
@constraint(model, ∂(x2) == u2)

obstacle = LazySets.Zonotope([-0.5, 0.0], [0.6 1.1; 1.3 -0.5])
target = LazySets.Ball2([3.0, 3.0], 0.8)

@constraint(model, start(x1) in MOI.Interval(-3.6, -3.4))
@constraint(model, start(x2) in MOI.Interval(-3.6, -3.4))
@constraint(model, [x1, x2] ∉ obstacle)
@constraint(model, [x1, x2] in Final(target))

set_attribute(model, "jacobian_bound", u -> SMatrix{2, 2}(0.0, 0.0, 0.0, 0.0))
set_attribute(model, "time_step", 0.3)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.2, 0.2)))
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))
set_attribute(model, "print_level", 0)
optimize!(model)
println("termination: ", termination_status(model))

trajectory = Dionysos.simulate(model, SVector(-3.5, -3.5); nsteps = 160)

concrete_problem = get_attribute(model, "concrete_problem")
fig = plot(; aspect_ratio = :equal, legend = false, dpi = 200)
plot!(fig, concrete_problem)
plot!(fig, obstacle; color = :black, opacity = 0.8)
plot!(fig, trajectory; ms = 1.2, color = :blue)
out = joinpath(@__DIR__, "zonotope_reach_avoid.png")
savefig(fig, out)
println("saved ", out)
