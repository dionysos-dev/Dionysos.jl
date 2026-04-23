using StaticArrays, LinearAlgebra, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

using MathematicalSystems
MS = MathematicalSystems

using JLD2

include("../../problems/Pendulum/simple_pendulum.jl")

system = SimplePendulum.system(;
    l = 1.0,
    g = 9.81,
    _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0)),
    _U_ = UT.HyperRectangle(SVector(-10.0), SVector(10.0)),
)

tstep = 0.1
discrete_time_system = ST.discretize_continuous_system(system, tstep)

# ------------------------------------------------------------
# PID Controller on Position
# ------------------------------------------------------------

input_set = MS.inputset(system)

pid_controller = ST.PIDControllers.PIDControllerVector(;
    Kp = SVector(15.0, 5.0),
    Ki = SVector(3.0, 0.0),
    Kd = SVector(0.0, 0.0),
    ref = ST.PIDControllers.ConstantSignal(SVector(3π / 4.0, 0.0)),
    error = ST.PIDControllers.WrapAnglePositionVelocityError(),
    dt = tstep,
    time_getter = ST.PIDControllers.ConstantTimeGetter(),
    umin = input_set.lb,
    umax = input_set.ub,
    e0 = SVector(0.0, 0.0),
)

# JLD2.jldsave("pid_controller.jld2"; controller = pid_controller)
# pid_controller = JLD2.load("pid_controller.jld2", "controller")

# ------------------------------------------------------------
# Closed-loop simulation
# ------------------------------------------------------------

nstep = 200
x0 = SVector(0.0, 0.0)

traj = ST.get_closed_loop_trajectory(discrete_time_system, pid_controller, x0, nstep)

x_traj = traj.x
u_traj = traj.u
# traj.q is the PID internal memory trajectory if your rollout returns it

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

fig = plot(; aspect_ratio = :equal)
plot!(system.X; color = :grey, hole_color = :black, opacity = 1.0, label = "")
plot!(x_traj; ms = 2.0, arrows = false)
display(fig)

# ------------------------------------------------------------
# Visualization
# ------------------------------------------------------------

using RigidBodyDynamics
using MeshCat, MeshCatMechanisms

urdf = joinpath(dirname(dirname(pathof(Dionysos))), "problems/pendulum/", "Pendulum.urdf")
mechanism = parse_urdf(urdf)
state = MechanismState(mechanism)
joint = first(joints(mechanism))

state_values = [x_traj.seq[i] for i in 1:ST.length(x_traj)]
ts = collect(0.0:tstep:((length(state_values) - 1) * tstep))
mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf))
vis = mvis.visualizer
open(vis)

fps = round(Int, 1 / tstep)
anim = MeshCat.Animation(vis; fps = fps)

for k in eachindex(ts)
    θ = state_values[k][1]
    set_configuration!(state, joint, θ)

    MeshCat.atframe(anim, k) do
        return MeshCatMechanisms.set_configuration!(mvis, configuration(state))
    end
end

MeshCat.setanimation!(vis, anim; play = true)
