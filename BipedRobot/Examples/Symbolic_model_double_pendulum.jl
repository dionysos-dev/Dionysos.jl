using Revise
using Latexify
using RigidBodyDynamics
using StaticArrays
using Symbolics

# ## Create symbolic parameters
# * Masses: $m_1, m_2$
# * Mass moments of inertia (about center of mass): $I_1, I_2$
# * Link lengths: $l_1, l_2$
# * Center of mass locations (w.r.t. preceding joint axis): $c_1, c_2$
# * Gravitational acceleration: $g$

inertias = @variables m_1 m_2 I_1 I_2 positive = true
lengths = @variables l_1 l_2 c_1 c_2 real = true
gravitational_acceleration = @variables g real = true
params = [inertias..., lengths..., gravitational_acceleration...]
transpose(params)

# ## Create double pendulum `Mechanism`
# A `Mechanism` contains the joint layout and inertia parameters, but no state information.

T = Num # the 'type' of the Mechanism we'll construct
axis = SVector(zero(T), one(T), zero(T)) # axis of rotation for each of the joints
double_pendulum = Mechanism(RigidBody{T}("world"); gravity = SVector(zero(T), zero(T), g))
world = root_body(double_pendulum) # the fixed 'world' rigid body

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

## The method for creating the QuaternionFloating joint was developed in my personal code
inertia1 = SpatialInertia(
    CartesianFrame3D("upper_link");
    moment = I_1 * axis * transpose(axis),
    com = SVector(zero(T), zero(T), c_1),
    mass = m_1,
)
body1 = RigidBody(inertia1)
joint1 = Joint("shoulder", Revolute(axis))
joint1_to_world = one(Transform3D{T}, frame_before(joint1), default_frame(world)); # basically this transform is the identity because the joint coordinate coincides with the origin of the world
attach!(double_pendulum, world, body1, joint1; joint_pose = joint1_to_world);

# Attach the second (lower) link to the world via a revolute joint named 'elbow'
inertia2 = SpatialInertia(
    CartesianFrame3D("lower_link");
    moment = I_2 * axis * transpose(axis),
    com = SVector(zero(T), zero(T), c_2),
    mass = m_2,
)
body2 = RigidBody(inertia2)
joint2 = Joint("elbow", Revolute(axis))
joint2_to_body1 =
    Transform3D(frame_before(joint2), default_frame(body1), SVector(zero(T), zero(T), l_1))
attach!(double_pendulum, body1, body2, joint2; joint_pose = joint2_to_body1)

## Export URDF
# write_urdf("test.urdf", double_pendulum; robot_name="double_pendulum", include_root=true)

# ## Create `MechanismState` associated with the double pendulum `Mechanism`
# A `MechanismState` stores all state-dependent information associated with a `Mechanism`.

x = MechanismState(double_pendulum);

# Set the joint configuration vector of the MechanismState to a new vector of symbolic variables

q = configuration(x)

q_s = Symbolics.variables(:q, 1:length(configuration(x)))
set_configuration!(x, q_s) # starting a pass initial configuration for i in eachindex(q)
q = configuration(x)

# Velocities vector
v = velocity(x)

v_s = Symbolics.variables(:v, 1:length(velocity(x)))
set_velocity!(x, v_s) # Set the joint velocity vector of the MechanismState to a new vector of symbolic variables
v = velocity(x)

# ## Compute dynamical quantities in symbolic form

# Mass matrix
MM = simplify.(mass_matrix(x)) # This gives you the general inertia matrix M(x) of the mechanism
# latexify(MM) #if we want it for latex

# Kinetic energy

simplify.(kinetic_energy(x))

# Potential energy

# this works in my dev version of RigidBodyDynamics because I eliminated the line with the boolean test in the function
simplify.(gravitational_potential_energy(x))

# Compute dynamics_bias c(q,v,w_ext)

simplify.(dynamics_bias(x))

# We can also do inverse dynamics...

v̇ = similar(v) # the joint acceleration vector, i.e., the time derivative of the joint velocity vector v
a_s = Symbolics.variables(:v̇, 1:length(velocity(x)))

for l in 1:length(v̇)
    v̇[l] = a_s[l]
end

simplify.(inverse_dynamics(x, v̇)) # this gives you τ

## M(q)v̇ +c(q,v,w) = τ, where q is the joint configuration vector (angles), v is the joint velocity vector, and v̇ is the joint. Furthermore, c(q,v,w) 
## is the Coriolis tensor, which embeds viscous friction torques and possible external signals (disturbances) w and effects of gravity. 
simplify.(inverse_dynamics(x, v̇)) + simplify.(dynamics_bias(x))
