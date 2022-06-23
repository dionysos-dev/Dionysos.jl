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

## modify the mechanism
# Attach the first (upper) link to the world via a revolute joint named 'shoulder'
inertia1 = SpatialInertia(CartesianFrame3D("upper_link"),
    moment=I_1 * axis * transpose(axis),
    com=SVector(zero(T), zero(T), c_1),
    mass=m_1)
body1 = RigidBody(inertia1)
joint1 = Joint("shoulder", Revolute(axis))
joint1_to_world = one(Transform3D{T}, frame_before(joint1), default_frame(world));
attach!(double_pendulum, world, body1, joint1,
    joint_pose = joint1_to_world);

# Attach the second (lower) link to the world via a revolute joint named 'elbow'
inertia2 = SpatialInertia(CartesianFrame3D("lower_link"),
    moment=I_2 * axis * transpose(axis),
    com=SVector(zero(T), zero(T), c_2),
    mass=m_2)
body2 = RigidBody(inertia2)
joint2 = Joint("elbow", Revolute(axis))
joint2_to_body1 = Transform3D(
    frame_before(joint2), default_frame(body1), SVector(zero(T), zero(T), l_1))
attach!(double_pendulum, body1, body2, joint2,
    joint_pose = joint2_to_body1)

# ## Create `MechanismState` associated with the double pendulum `Mechanism`
# A `MechanismState` stores all state-dependent information associated with a `Mechanism`.

x = MechanismState(double_pendulum);

# Joints vector
q = configuration(x)

# Velocities vector
v = velocity(x)

# ## Compute dynamical quantities in symbolic form

# Mass matrix
simplify.(mass_matrix(x)) # This gives you the general inertia matrix M(x) of the mechanism

# We can also do inverse dynamics...
@variables v_1 v_2 # joint velocities
@variables v̇_1 v̇_2 # joint accelarations

v̇ = similar(velocity(x)) # the joint acceleration vector, i.e., the time derivative of the joint velocity vector v
v̇[joint1][1] = v̇_1
v̇[joint2][1] = v̇_2

simplify.(inverse_dynamics(x, v̇)) # this gives you τ

## M(x)v̇ +c(q,v,w) = τ, where x is the joint configuration vector (angles), v is the joint velocity vector, and v̇ is the joint. Furthermore, c(q,v,w) 
## is the Coriolis tensor, which embeds viscous friction torques and possible external signals (disturbances) w. In this example, there is no c(q,v,w). 
