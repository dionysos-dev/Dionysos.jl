# Here we hope to chieve simulation of our biped robot. Options with and without control added.
# using Revise
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using Random
using StaticArrays
using Symbolics
using BipedRobot

## first we will try do define the robot mechanism with symbolic variables. Then we will export the robot model. 
## there are some limitations with this method. for example, only inertia tags are supported for the links
## second we will replace the symbolic variables by float numbers and export the URDF
## third we will add contact points and simulate the robot  

## First step
# inertias = @variables m_hip m_r_thigh m_r_leg m_r_foot m_l_thigh m_l_leg m_l_foot  I_hip I_r_thigh I_r_leg I_r_foot I_l_thigh I_l_leg I_l_foot positive = true
inertias =
    @variables m_h m_r_t m_r_l m_r_f m_l_t m_l_l m_l_f I_h I_r_t I_r_l I_r_f I_l_t I_l_l I_l_f positive =
        true
lengths =
    @variables l_h_x l_h_y l_h_z l_r_t c_r_t l_r_l c_r_l c_r_f l_l_t c_l_t l_l_l c_l_l c_l_f real =
        true
## the variables are, in order, the x,y,z lengths of the hips (a rectangular box for the moment), 
## the lenght and circunference of right thigh and of the right leg (cilindrical for the moment)
## the lenght and circunference of left thigh and of the left leg (cilindrical for the moment)

gravitational_acceleration = @variables g real = true
params = [inertias..., lengths..., gravitational_acceleration...]
transpose(params)

# ## Create robot `Mechanism`
# A `Mechanism` contains the joint layout and inertia parameters, but no state information.

T = Num # the 'type' of the Mechanism we'll construct
axis = SVector(zero(T), one(T), zero(T)) # axis of rotation for each for the knee and hips joints (around y-axis)
robot = Mechanism(RigidBody{T}("world"); gravity = SVector(zero(T), zero(T), g))
world = root_body(robot) # the fixed 'world' rigid body

## modify the mechanism

# create a symbolic symmetric inertia matrix
function get_symb_Inertia_matrix(H)
    I_m = Symbolics.variables(H, 1:3, 1:3)
    for i in 1:3
        for j in 1:3
            if i > j
                I_m[i, j] = I_m[j, i]
            end
        end
    end
    return I_m
end

"""

Construct a `SpatialInertia` by specifying:

* `frame`: the frame in which the spatial inertia is expressed.
* one of:
  * `moment`: the moment of inertia expressed in `frame` (i.e., about the origin of `frame` and in `frame`'s axes).
  * `moment_about_com`: the moment of inertia about the center of mass, in `frame`'s axes.
* `com`: the center of mass expressed in `frame`.
* `mass`: the total mass.

The `com` and `mass` keyword arguments are required, as well as exactly one of `moment` and `moment_about_com`
"""

## Attach hips

inertia1 = SpatialInertia(
    CartesianFrame3D("hips_link");
    moment = Symbolics.variables(:I_h, 1:3, 1:3),
    com = SVector(l_h_x, l_h_y, l_h_z),
    mass = m_h,
)

body1 = RigidBody(inertia1)

# Attach the hip link to the world via a floating joint named 'world_to_hips'

joint1 = Joint("world_to_hips", QuaternionFloating{Num}())

# joint1_to_world = one(Transform3D{T}, frame_before(joint1), default_frame(world));
joint1_to_world = Transform3D(
    frame_before(joint1),
    default_frame(world),
    SVector(zero(T), one(T), one(T)),
)

attach!(robot, world, body1, joint1; joint_pose = joint1_to_world); # attach the link to the mechanism tree

# We now repeat the process to the other parts of the robot.

## Attach the right thigh link to the hips link via a revolute joint named 'r_hips_joint'
inertia2 = SpatialInertia(
    CartesianFrame3D("r_thigh_link");
    moment = I_r_t * axis * transpose(axis),
    com = SVector(zero(T), zero(T), c_r_t),
    mass = m_r_t,
)

body2 = RigidBody(inertia2)

joint2 = Joint("r_hips_joint", Revolute(axis))

joint2_to_body1 = Transform3D(
    frame_before(joint2),
    default_frame(body1),
    SVector(zero(T), zero(T), -l_r_t),
)

attach!(robot, body1, body2, joint2; joint_pose = joint2_to_body1)

## Attach the left thigh link to the hips link via a revolute joint named 'l_hips_joint'
inertia4 = SpatialInertia(
    CartesianFrame3D("l_thigh_link");
    moment = I_l_t * axis * transpose(axis),
    com = SVector(zero(T), zero(T), c_l_t),
    mass = m_l_t,
)

body4 = RigidBody(inertia4)

joint4 = Joint("l_hips_joint", Revolute(axis))

joint4_to_body1 = Transform3D(
    frame_before(joint4),
    default_frame(body1),
    SVector(zero(T), zero(T), -l_l_t),
)

attach!(robot, body1, body4, joint4; joint_pose = joint4_to_body1)

## Attach the right leg link to the right thigh link via a revolute joint named 'r_knee_joint'
inertia3 = SpatialInertia(
    CartesianFrame3D("r_leg_link");
    moment = I_r_l * axis * transpose(axis),
    com = SVector(zero(T), zero(T), c_r_l),
    mass = m_r_l,
)

body3 = RigidBody(inertia3)

joint3 = Joint("r_knee_joint", Revolute(axis))

joint3_to_body2 = Transform3D(
    frame_before(joint3),
    default_frame(body2),
    SVector(zero(T), zero(T), -l_r_l),
)

attach!(robot, body2, body3, joint3; joint_pose = joint3_to_body2)

## Attach the left leg link to the left thigh link via a revolute joint named 'l_knee_joint'
inertia5 = SpatialInertia(
    CartesianFrame3D("l_leg_link");
    moment = I_l_l * axis * transpose(axis),
    com = SVector(zero(T), zero(T), c_l_l),
    mass = m_l_l,
)

body5 = RigidBody(inertia5)

joint5 = Joint("l_knee_joint", Revolute(axis))

joint5_to_body4 = Transform3D(
    frame_before(joint5),
    default_frame(body4),
    SVector(zero(T), zero(T), -l_r_l),
)

attach!(robot, body4, body5, joint5; joint_pose = joint5_to_body4)

# ## to be continued...

## Export URDF
#write_urdf("biped_robot_exported.urdf", robot; robot_name="biped_robot_exported", include_root=true)
## the problem with the write_urdf function is that is does not include some tags as the visual one

# # ## Create `MechanismState` associated with the double pendulum `Mechanism`
# # A `MechanismState` stores all state-dependent information associated with a `Mechanism`.

# x = MechanismState(robot);

# # Set the joint configuration vector of the MechanismState to a new vector of symbolic variables

# q = configuration(x)

# q_s = Symbolics.variables(:q, 1:length(configuration(x)))
# set_configuration!(x, q_s) # starting a pass initial configuration for i in eachindex(q)
# q = configuration(x)

# # Velocities vector
# v = velocity(x)

# v_s = Symbolics.variables(:v, 1:length(velocity(x)))
# set_velocity!(x, v_s) # Set the joint velocity vector of the MechanismState to a new vector of symbolic variables
# v = velocity(x)

# # ## Compute dynamical quantities in symbolic form

# # Mass matrix
# MM = simplify.(mass_matrix(x)) # This gives you the general inertia matrix M(x) of the mechanism
# # latexify(MM) #if we want it for latex

# # Kinetic energy

# simplify.(kinetic_energy(x))

# # Potential energy

# # this works in my dev version of RigidBodyDynamics because I eliminated the line with the boolean test in the function
# simplify.(gravitational_potential_energy(x))

# # Compute dynamics_bias c(q,v,w_ext)

# simplify.(dynamics_bias(x))

# # We can also do inverse dynamics...

# v̇ = similar(v) # the joint acceleration vector, i.e., the time derivative of the joint velocity vector v
# a_s = Symbolics.variables(:v̇, 1:length(velocity(x)))

# for l = 1:length(v̇)
#     v̇[l] = a_s[l]
# end

# simplify.(inverse_dynamics(x, v̇)) # this gives you τ

# ## M(q)v̇ +c(q,v,w) = τ, where q is the joint configuration vector (angles), v is the joint velocity vector, and v̇ is the joint. Furthermore, c(q,v,w) 
# ## is the Coriolis tensor, which embeds viscous friction torques and possible external signals (disturbances) w and effects of gravity. 
# simplify.(inverse_dynamics(x, v̇))+simplify.(dynamics_bias(x))

# #### If we want to load the URDF instead of creating the robot with RigidBodyDynamics
## describe path to urdf file
packagepath() = joinpath(@__DIR__, "..", "deps")
urdfpath() = joinpath(packagepath(), "biped_robot.urdf")

## load mechanism 'robot' from the URDF file
robot = parse_urdf(Float64, urdfpath())

# list_of_joints = joints(robot)
# list_of_links = bodies(robot)

## here create function to read all the bodies and change their masses, center of mass, and moment of inertia to symbolic variables in order to obtain model
# for body in bodies(robot)
#     l = 0
#     if body.name != "world"
#     body.inertia.mass = inertias[l+1] # mass of the link
#     body.inertia.cross_part = SVector(l_h_x,l_h_y,l_h_z) # center of mass of the link
#     body.inertia.moment = I_h * axis * transpose(axis); #moment of mass matrix (3 by 3 SMatrix)
#     l = l+1
#     end
# end

## add flat ground or not
add_flat_ground = true

if add_flat_ground
    frame = root_frame(robot)
    ground = HalfSpace3D(Point3D(frame, 0.0, 0.0, 0.0), FreeVector3D(frame, 0.0, 0.0, 1.0))
    add_environment_primitive!(robot, ground)
end

## add contact points in the links

# define contact model
function default_contact_model()
    return SoftContactModel(
        hunt_crossley_hertz(; k = 500e3),
        ViscoelasticCoulombModel(0.8, 20e3, 100.0),
    )
end
contactmodel = default_contact_model()

# # create contact point in right leg_link
right_foot_link = findbody(robot, "r_foot_link")
frame_lower_link = default_frame(right_foot_link)
add_contact_point!(
    right_foot_link,
    ContactPoint(Point3D(frame_lower_link, 0.0, 0.0, -0.03), contactmodel),
)

# # create contact point in left leg_link
left_foot_link = findbody(robot, "l_foot_link")
frame_lower_link = default_frame(left_foot_link)
add_contact_point!(
    left_foot_link,
    ContactPoint(Point3D(frame_lower_link, 0.0, 0.0, -0.03), contactmodel),
)

## the state of the mechanism is a type that contains many information about the mechanism. most importantly, the angles, velocities, etc of the joints
const state = MechanismState(robot)

## visualisation and simulation

@static if get(ENV, "CI", "false") == "false"
    using MeshCatMechanisms
    vis = MechanismVisualizer(robot, URDFVisuals(urdfpath()))
    open(vis)
end

## set the configurations and velocities of the joints (i.e., initial angles (called configuration in julia robotics) and initial velocities):
# set_configuration!(state, [1,0,0,0,0,0,1,0,0,0,0]) # starting a pass initial configuration
set_configuration!(state, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]) # starting a pass initial configuration

@static if get(ENV, "CI", "false") == "false"
    set_configuration!(vis, configuration(state)) ## update the configuration also in the visualiser
end

# ## Basic simulation is easy (but see RigidBodySim.jl for a more featureful simulator). 

# ## This simulation option has no control
# ts, qs, vs = simulate(state, 0.22, Δt = 1e-3);

# ## with control. This serves to illustrate how to create feedback controllers in the simulation

controllable_joints = 7:10

integrator = zeros(10)

function control!(torques::AbstractVector, t, state::MechanismState)
    # rand!(torques) # for example
    for joint in joints(robot)
        v_range = velocity_range(state, joint)
        # print(l)
        # torques[velocity_range(state, joint)] .= -1.0*velocity(state, joint) # feedbacking the velocity in each joint
        # if configuration_range(state, joint) < 11:11  # feedbacking the angle in each joint. using the if because there is one more angle than torques in this atls robot, which I don't know why
        # torques[configuration_range(state, joint)] .= -1*configuration(state,joint)
        # end
        if !isempty(v_range) && v_range ⊆ controllable_joints
            print(state.q[v_range])
            print("\t")
            torques[v_range] .=
                -50 * state.q[v_range] - 1 * state.v[v_range] + 10 * integrator[v_range]
            integrator[v_range] .= integrator[v_range] - 0.001 * state.q[v_range]
        else
            torques[v_range] .= 0 # no control
        end
    end
    return println()
end

ts, qs, vs = simulate(state, 1, control!; Δt = 1e-3);

## After which we can animate the results:
@static if get(ENV, "CI", "false") == "false"
    MeshCatMechanisms.animate(vis, ts, qs; realtimerate = 0.2)
end

##########
