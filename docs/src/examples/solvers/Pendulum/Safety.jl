using Test
using StaticArrays, Plots

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(joinpath(dirname(dirname(pathof(Dionysos))), "problems/pendulum", "pendulum.jl"))

concrete_problem = Pendulum.problem(; objective = "safety-up")
concrete_system = concrete_problem.system

x0 = SVector(0.0, 0.0)
hx = SVector(0.05, 0.05)
state_grid = DO.GridFree(x0, hx)

u0 = SVector(0.0);
h = SVector(0.3);
input_grid = DO.GridFree(u0, h);

using JuMP
optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.optimize!(optimizer)

abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

# ### Trajectory display
nstep = 100
x0 = SVector(UT.sample(concrete_problem.initial_set)...)

control_trajectory =
    ST.get_closed_loop_trajectory(concrete_system.f, concrete_controller, x0, nstep)

# fig = plot(; aspect_ratio = :equal);
# plot!(concrete_system.X);
# plot!(concrete_problem.safe_set; markersize = 1, arrows = false)
# plot!(control_trajectory; markersize = 1, arrows = false)
# display(fig)

fig = plot(; aspect_ratio = :equal, legend=false);
light_grey = "#E8E8E8"
grey_color = "#CCCCCC"
_O_ = UT.HyperRectangle(SVector(-2.9, -10.0), SVector(-2.1, 10.0)) #2.9
plot!(concrete_system.X;color=light_grey);
# plot!(abstract_system.Xdom; color = grey_color, opacity = 1.0);
plot!(concrete_problem.safe_set; dims = [1, 2], color = :red, opacity = 0.8);
plot!(concrete_problem.initial_set; dims = [1, 2], color = :green, opacity = 0.8);
# plot!(concrete_problem.target_set; dims = [1, 2], color = :red, opacity = 0.8);
plot!(_O_; dims = [1, 2], color = :black, opacity = 0.8);
plot!(control_trajectory; color=:blue, markersize = 1, arrows = false)
xlims!(fig, (-pi-0.5 , pi+0.5))
ylims!(fig, (-10.5, 10.5))
display(fig)

# savefig("pendulum_safety_problem.pdf")
savefig("pendulum_safety_traj.pdf")



# ### For Visualization
using RigidBodyDynamics, MeshCatMechanisms
urdf = joinpath(dirname(dirname(pathof(Dionysos))), "problems/pendulum/", "Pendulum.urdf")
mechanism = parse_urdf(urdf)
state = MechanismState(mechanism)
joint = first(joints(mechanism))

tstep = 0.1
state_values =
    [ST.get_state(control_trajectory, i)[1] for i in 1:ST.length(control_trajectory)]
ts = collect(0.0:tstep:((length(state_values) - 1) * tstep))
qs = Vector{Vector{Float64}}(undef, length(state_values))
for i in 1:length(state_values)
    set_configuration!(state, joint, state_values[i])
    qs[i] = configuration(state)
end
mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf))
# Define the obstacle as a box (HyperRectangle) and add it to the viewer
box_width, box_height, box_depth = 0.9, 0.5, 0.5
obstacle_geometry = GeometryBasics.HyperRectangle(Vec(0.0, 0.0, 0.0), Vec(box_width, box_height, box_depth))
setobject!(mvis["obstacle"], obstacle_geometry, MeshCat.MeshPhongMaterial(color=RGB(0.0, 0.0, 0.0)))
settransform!(mvis["obstacle"], Translation(0.2, -0.25, 0.43)) # left,0,h

MeshCatMechanisms.setanimation!(mvis, ts, qs)
