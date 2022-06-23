using RigidBodyDynamics
using RigidBodyDynamics.Contact
using Random
using StaticArrays
using MathematicalSystems
using Dionysos

# #### If we want to load the URDF instead of creating the robot with RigidBodyDynamics
## describe path to urdf file
#packagepath() = joinpath(@__DIR__, "..","..", "deps")
urdfpath() = joinpath(@__DIR__, "biped_robot.urdf")

## load mechanism 'robot' from the URDF file
robot = parse_urdf(Float64, urdfpath())

## add flat ground or not
add_flat_ground=true

if add_flat_ground
    frame = root_frame(robot)
    ground = HalfSpace3D(Point3D(frame, 0., 0., 0.), FreeVector3D(frame, 0., 0., 1.))
    add_environment_primitive!(robot, ground)
end

## add contact points in the links

# define contact model
function default_contact_model()
    SoftContactModel(hunt_crossley_hertz(k = 500e3), ViscoelasticCoulombModel(0.8, 20e3, 100.))
end
contactmodel = default_contact_model()

# # create contact point in right leg_link
right_leg_link = findbody(robot,"r_leg_link")
frame_lower_link = default_frame(right_leg_link)
add_contact_point!(right_leg_link, ContactPoint(Point3D(frame_lower_link, 0.0, 0.0, -1.05), contactmodel))

# # create contact point in left leg_link
left_leg_link = findbody(robot,"l_leg_link")
frame_lower_link = default_frame(left_leg_link)
add_contact_point!(left_leg_link, ContactPoint(Point3D(frame_lower_link, 0.0, 0.0, -1.05), contactmodel))


## the state of the mechanism is a type that contains many information about the mechanism. most importantly, the angles, velocities, etc of the joints
const state = MechanismState(robot)

set_configuration!(state, [1,0,0,0,0,0,0,0,0,0,-pi/4]) # starting a pass initial configuration

function vectorFieldBipedRobot(x,u)
    q=x[1:num_positions(state)]
    v=x[num_positions(state)+1:num_positions(state)+num_velocities(state)]
    s=x[num_positions(state)+num_velocities(state)+1:end]
    
    set_configuration!(state,q)
    set_velocity!(state,v)
    set_additional_state!(state,s)
    
    ẋ=similar(x)
    
    dynamics!(ẋ,DynamicsResult(state.mechanism),state,x,u)
    
    ẋ
end

system = BlackBoxControlContinuousSystem(vectorFieldBipedRobot,
    num_positions(state)+num_velocities(state)+num_additional_states(state),
    num_velocities(state))

# TODO: create the Dionysos.Control.OptimalControlProblem