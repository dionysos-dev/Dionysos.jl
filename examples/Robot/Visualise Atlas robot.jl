# This example shows how to load and simulate the Atlas robot. No feedback control at this moment.

using AtlasRobot
using RigidBodyDynamics
using MeshCatMechanisms
using RigidBodyDynamics.Contact: location

# mechanism = AtlasRobot.mechanism(add_flat_ground=true,contactmodel=nothing) # in this case there is no contact and robot falls in the world
mechanism = AtlasRobot.mechanism(add_flat_ground=true)

    # mechanism = AtlasRobot.mechanism()
    # meshdir = joinpath(AtlasRobot.packagepath(), "Atlas", "urdf", "meshes")
    state = MechanismState(mechanism)
    AtlasRobot.setnominal!(state)
    for sideprefix in ('l', 'r')
        foot = findbody(mechanism, "$(sideprefix)_foot")
        for point in contact_points(foot)
            contact_location_world = transform(state, location(point), root_frame(mechanism))
        end
    end


## visualisation
vis = MechanismVisualizer(mechanism);
open(vis)

## set the configurations and velocities of the joints:
#set_configuration!(state, [1.0, -1.5])
set_configuration!(vis, configuration(state))


## Basic simulation is easy (but see RigidBodySim.jl for a more featureful simulator):

ts, qs, vs = simulate(state, 5., Δt = 1e-3);

## After which we can animate the results:
MeshCatMechanisms.animate(vis, ts, qs; realtimerate = 0.25)

##########
