using StaticArrays, JuMP, LinearAlgebra, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

using MathematicalSystems
MS = MathematicalSystems

include("../problems/simple_pendulum.jl");

system = SimplePendulum.system(;
    l = 1.0,
    g = 9.81,
    _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0)),
    _U_ = UT.HyperRectangle(SVector(-10), SVector(10)),
)

tstep = 0.1
wrap_angle(e) = mod(e + π, 2π) - π
error = (x, r, _) -> begin
    eθ = wrap_angle(r[1] - x[1])
    eω = r[2] - x[2]
    return SVector(eθ, eω)
end

input_set = MS.inputset(system)
pid = ST.PIDControllers.PIDControllerVector(;
    Kp = SVector(15.0, 5.0), # SVector(20.0, 5.0)
    Ki = SVector(3.0, 0.0),  # SVector(1.0, 0.0)
    Kd = SVector(0.0, 0.0),
    ref = SVector(3*pi/4.0, 0.0),
    error = error,
    dt = tstep,
    umin = input_set.lb,
    umax = input_set.ub,
    e0 = SVector(0.0, 0.0),
)
pid_controller = ST.PIDControllers.pid_map(pid; nx = 2, nu = 1, silent = false)

discrete_time_system = ST.discretize_continuous_system(system, tstep)

nstep = 200
x0 = SVector(0.0, 0.0)
x_traj, u_traj =
    ST.get_closed_loop_trajectory(discrete_time_system, pid_controller, x0, nstep;)

fig = plot(; aspect_ratio = :equal);
plot!(system.X; color = :grey, hole_color = :black, opacity = 1.0, label = "");
plot!(x_traj; ms = 2.0, arrows = false)
display(fig)

# ### For Visualization

using RigidBodyDynamics
using MeshCat, MeshCatMechanisms

# --- build mechanism/state from your URDF ---
urdf = joinpath(dirname(dirname(pathof(Dionysos))), "problems/pendulum/", "Pendulum.urdf")
mechanism = parse_urdf(urdf)
state = MechanismState(mechanism)
joint = first(joints(mechanism))

# --- build trajectory data ---
state_values = [x_traj.seq[i] for i in 1:ST.length(x_traj)]
ts = collect(0.0:tstep:((length(state_values) - 1) * tstep))

# --- visualizer ---
mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf))
vis = mvis.visualizer
open(vis)

# --- animation (no Interpolations) ---
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
