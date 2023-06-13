using RigidBodyDynamics
using StaticArrays
using Symbolics

# ## Create symbolic parameters
# * Masses: $m_1, m_2$
# * Mass moments of inertia (about center of mass): $I_1, I_2$
# * Link lengths: $l_1, l_2$
# * Center of mass locations (w.r.t. preceding joint axis): $c_1, c_2$
# * Gravitational acceleration: $g$

inertias = @variables m_1 I_1 positive = true
lengths = @variables l_1 c_1 real = true
gravitational_acceleration = @variables g real = true
params = [inertias..., lengths..., gravitational_acceleration...]
transpose(params)

# ## Create double pendulum `Mechanism`
# A `Mechanism` contains the joint layout and inertia parameters, but no state information.

T = Num # the 'type' of the Mechanism we'll construct
axis = SVector(zero(T), one(T), zero(T)) # axis of rotation for each of the joints
double_pendulum = Mechanism(RigidBody{T}("world"); gravity = SVector(zero(T), zero(T), g))
world = root_body(double_pendulum) # the fixed 'world' rigid body

## modify the mechanism
# Attach the first (upper) link to the world via a revolute joint named 'shoulder'
inertia1 = SpatialInertia(
    CartesianFrame3D("upper_link");
    moment = I_1 * axis * transpose(axis),
    com = SVector(zero(T), zero(T), c_1),
    mass = m_1,
)
body1 = RigidBody(inertia1)
joint1 = Joint("shoulder", Revolute(axis))
joint1_to_world = one(Transform3D{T}, frame_before(joint1), default_frame(world));
attach!(double_pendulum, world, body1, joint1; joint_pose = joint1_to_world);

## Export URDF
# write_urdf("test.urdf", double_pendulum; robot_name="double_pendulum", include_root=true)

# ## Create `MechanismState` associated with the double pendulum `Mechanism`
# A `MechanismState` stores all state-dependent information associated with a `Mechanism`.

x = MechanismState(double_pendulum);

# Set the joint configuration vector of the MechanismState to a new vector of symbolic variables

q = configuration(x)

q_s = @variables q_1
set_configuration!(x, q_s) # starting a pass initial configuration for i in eachindex(q)
q = configuration(x)

# Velocities vector
v = velocity(x)
v_s = @variables v_1
set_velocity!(x, v_s) # Set the joint velocity vector of the MechanismState to a new vector of symbolic variables
v = velocity(x)

# ## Compute dynamical quantities in symbolic form

# Mass matrix
simplify.(mass_matrix(x)) # This gives you the general inertia matrix M(x) of the mechanism

# Kinetic energy

simplify.(kinetic_energy(x))

# Potential energy

# this works in my dev version of RigidBodyDynamics because I eliminated the line with the boolean test in the function
simplify.(gravitational_potential_energy(x))

# Compute dynamics_bias c(q,v,w_ext)

simplify.(dynamics_bias(x))

# We can also do inverse dynamics...
@variables v̇_1 # joint accelarations

v̇ = similar(v) # the joint acceleration vector, i.e., the time derivative of the joint velocity vector v
v̇[joint1][1] = v̇_1

simplify.(inverse_dynamics(x, v̇)) # this gives you τ

## M(x)v̇ +c(q,v,w) = τ, where x is the joint configuration vector (angles), v is the joint velocity vector, and v̇ is the joint. Furthermore, c(q,v,w) 
## is the Coriolis tensor, which embeds viscous friction torques and possible external signals (disturbances) w and effects of gravity. 
simplify.(inverse_dynamics(x, v̇)) + simplify.(dynamics_bias(x))
