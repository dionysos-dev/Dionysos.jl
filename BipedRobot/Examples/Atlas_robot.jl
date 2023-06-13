# This example shows how to load and simulate the Atlas robot. Options with and without control added.

using AtlasRobot
using RigidBodyDynamics
using MeshCatMechanisms
using RigidBodyDynamics.Contact: location
using Random

# mechanism = AtlasRobot.mechanism(add_flat_ground=true,contactmodel=nothing) # in this case there is no contact and robot falls in the world
mechanism = AtlasRobot.mechanism(; add_flat_ground = true)
# mechanism = AtlasRobot.mechanism()  # in this case there is also no contact because the ground is not added and robot falls in the world

state = MechanismState(mechanism)
AtlasRobot.setnominal!(state)

## test export urdf
#write_urdf("test.urdf", mechanism; robot_name="atlas_robot", include_root=true)
## the problem with the write_urdf function is that is does not include some tags as the visual one

## visualisation

vis = MechanismVisualizer(mechanism);
open(vis)

## set the configurations and velocities of the joints:

set_configuration!(vis, configuration(state))

## Basic simulation is easy (but see RigidBodySim.jl for a more featureful simulator). 

## This simulation option has no control
# ts, qs, vs = simulate(state, 0.22, Δt = 1e-3);

## With control. This serves to illustrate how to create feedback controllers in the simulation.  
## The controllers however do not stabilise the Atlas Robot. 

function control!(torques::AbstractVector, t, state::MechanismState)
    # rand!(torques) # for example
    for joint in joints(mechanism)
        torques[velocity_range(state, joint)] .= 0 * velocity(state, joint) # feedbacking the velocity in each joint
        # if configuration_range(state, joint) < 37:37  # feedbacking the angle in each joint. using the if because there is one more angle than torques in this atls robot, which I don't know why
        # torques[configuration_range(state, joint)] .= -0.1*configuration(state,joint)
        # end
        # torques[velocity_range(state, joint)] .= 0 # no control
    end
end

ts, qs, vs = simulate(state, 1, control!; Δt = 1e-3);

## After which we can animate the results:
MeshCatMechanisms.animate(vis, ts, qs; realtimerate = 0.1)

##########
