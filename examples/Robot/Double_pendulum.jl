using RigidBodyDynamics
using Plots
using MeshCatMechanisms


# urdf = joinpath(Dyonisos.packagepath(), "/deps/BipedRobot", "/doublependulum.urdf")


packagepath() = joinpath(@__DIR__, "..","..", "deps")
urdfpath() = joinpath(packagepath(), "BipedRobot", "doublependulum.urdf")
# doublependulum = parse_urdf(Float64, urdf)


# urdf = "doublependulum.urdf"
doublependulum = parse_urdf(Float64, urdfpath())

const state = MechanismState(doublependulum)

vis = MechanismVisualizer(doublependulum, URDFVisuals(urdfpath()));
open(vis)

## set the configurations and velocities of the joints:
set_configuration!(state, [1.0, -1.5])
set_configuration!(vis, configuration(state))

## Basic simulation is easy (but see RigidBodySim.jl for a more featureful simulator):

ts, qs, vs = simulate(state, 5., Δt = 1e-3);

## After which we can animate the results:
MeshCatMechanisms.animate(vis, ts, qs)

## Or plot them using e.g. Plots.jl:
shoulder_angles = collect(q[1] for q in qs)
plot(ts, shoulder_angles, 
    xlabel = "Time [s]", 
    ylabel = "Angle [rad]", 
    label = "Shoulder")
