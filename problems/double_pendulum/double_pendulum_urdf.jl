module DoublePendulum
using RigidBodyDynamics, MeshCat, MeshCatMechanisms 
using LinearAlgebra, StaticArrays, Random, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

############### to generate the URDF file of a double pendulum ###############
# The Mechanism object represents the kinematic and dynamic properties of the robot model.

# g : gravitational acceleration in z-direction
# I1 : moment of inertia about joint axis
# lc1 : center of mass location with respect to joint axis
# m1 : mass
# l1 : length of body 1
# I2 : moment of inertia about joint axis
# lc2 : center of mass location with respect to joint axis
# m2 : mass
# l2 : length of body 2. The mechanism does not store the ‘length’ of the second link (since that’s not necessary for the dynamics),
# but for the printing
# impossible to write <visual> and <collision> tags due to the fact that Mechanisms do not store the required information to write these tags.
function create_pendulum_mechanism(g, I1, lc1, m1, l1, I2, lc2, m2)
    # add first body (world), creating a 'root' rigid body, representing the fixed world
    base_link = RigidBody{Float64}("base_link")
    mechanism = Mechanism(base_link; gravity = SVector(0, 0, -g))

    # add second body (first arm)
    axis = SVector(0., 1., 0.) # joint axis
    inertia1 = SpatialInertia(CartesianFrame3D("link1"), moment=I1 * axis * axis', com=SVector(0, 0, -lc1), mass=m1)
    link1 = RigidBody(inertia1)
    joint1 = Joint("joint1", Revolute(axis))
    joint1_to_base_link = one(Transform3D, joint1.frame_before, default_frame(base_link))
    attach!(mechanism, base_link, link1, joint1, joint_pose = joint1_to_base_link)

    # add third body (second arm)
    inertia2 = SpatialInertia(CartesianFrame3D("link2"), moment=I2 * axis * axis', com=SVector(0, 0, -lc2), mass=m2)
    link2 = RigidBody(inertia2)
    joint2 = Joint("joint2", Revolute(axis))
    joint2_to_link1 = Transform3D(joint2.frame_before, default_frame(link1), SVector(0, 0, -l1))
    attach!(mechanism, link1, link2, joint2, joint_pose = joint2_to_link1)

    return mechanism
end

function add_visual!(urdf_file, lc1, l1, l2; r = 0.05, bs = [0.4 0.3 0.1])
    # Define the visual properties of the bodies
    # 1: base_link
    new_link_tag = """
        <link name="base_link">
         <visual>
            <geometry>
            <box size="$(bs[1]) $(bs[2]) $(bs[3])" />
            </geometry>
            <material name="green">
            <color rgba="0 1 0 1" />
            </material>
        </visual>
        <collision name="base_link_collision">
            <geometry>
            <box size="$(bs[1]) $(bs[2]) $(bs[3])" />
            </geometry>
        </collision>
        </link>
        """
    UT.replace_occurrence!(urdf_file, "<link name=\"base_link\"/>", new_link_tag, 1)
    # 2: link1
    new_str1 = """
    <visual>
    <origin xyz="0 0 $(-lc1)" rpy="0 0 0" />
    <geometry>
        <cylinder length="$(l1)" radius="$r" />
    </geometry>
    <material name="red">
        <color rgba="1 0 0 1" />
    </material>
    </visual>
    <collision name="link1_collision">
    <origin xyz="0 0 $(-lc1)" rpy="0 0 0" />
    <geometry>
        <cylinder length="$(l1)" radius="$r" />
    </geometry>
    </collision>
    """
    UT.write_after_string!(urdf_file, "</inertial>", new_str1, 1)
    # 3: link2
    new_str2 = """
    <visual>
    <origin xyz="0 0 $(-l1)" rpy="0 0 0" />
    <geometry>
        <cylinder length="$(l2)" radius="$r" />
    </geometry>
    <material name="blue">
        <color rgba="0 0 1 1" />
    </material>
    </visual>
    <collision name="link2_collision">
    <origin xyz="0 0 $(-l1)" rpy="0 0 0" />
    <geometry>
        <cylinder length="$(l2)" radius="$r" />
    </geometry>
    </collision>
    """
    UT.write_after_string!(urdf_file, "</inertial>", new_str2, 2)
end
    
function create_urdf_file(urdf_file; g=9.81, I1=0.333, lc1=0.5, m1=1., l1=1., I2=1.33, lc2=1., m2=1., l2=2.)
    mechanism = create_pendulum_mechanism(g, I1, lc1, m1, l1, I2, lc2, m2)
    write_urdf(urdf_file, mechanism)
    add_visual!(urdf_file, lc1, l1, l2)
end

function set_configuration(mechanism, q, v)
    # Get the states
    states = MechanismState(mechanism)
    # Get the joints 
    jointsL = joints(mechanism)

    # set the configurations and velocities of the joints
    RigidBodyDynamics.set_configuration!(states, jointsL[1], q[1])
    RigidBodyDynamics.set_configuration!(states, jointsL[2], q[2])
    set_velocity!(states, jointsL[1], v[1])
    set_velocity!(states, jointsL[2], v[2])

    # **Important**: a `MechanismState` contains cache variables that depend on the configurations and velocities of the joints. 
    # These need to be invalidated when the configurations and velocities are changed. To do this, call
    setdirty!(states)
    return states
end

#control! will be called at each time step of the simulation and allows you to specify joint torques
#given the time and the state of the Mechanism.
#torques is a vector of size equal to the number of joint.
function control!(torques::AbstractVector, t, state::MechanismState)
    torques[1] = 0.0
    torques[2] = 0.0
    #rand!(torques) # for example
    #println(torques)
end

function simulate_double_pendulum(urdf)
    mechanism = parse_urdf(urdf)
    # simulation
    q = [π,Float64(0.1)]  #[π,Float64(0.1)]
    v = [0.0,0.0]
    states = MechanismState(mechanism, q, v)
    #states = set_configuration(mechanism, q, v)
    ts, qs, vs = simulate(states, 5.0, control!; Δt = 1e-3)

    # visualization
    vis = Visualizer()
    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf), vis)

    animation = MeshCat.Animation(mvis, ts, qs)
    setanimation!(mvis, animation)
    open(vis)  # open the visualizer in a separate tab/window
end


function plot_double_pendulum(urdf)
    mechanism = parse_urdf(urdf)
    # simulation
    q = [π,Float64(0.1)]
    v = [0.0,0.0]
    states = MechanismState(mechanism, q, v)
    #states = set_configuration(mechanism, q, v)
    ts, qs, vs = simulate(states, 5.0, control!; Δt = 1e-3)

    p = plot(xlabel = "Time [s]", ylabel = "Angle [rad]")
    plot!(ts, collect(q[1] for q in qs), label = "Shoulder") 
    plot!(ts, collect(q[2] for q in qs), label = "Elbow") 
    plot!(ts, collect(v[1] for v in vs), label = "Derivative Shoulder") 
    plot!(ts, collect(v[2] for v in vs), label = "Derivative Elbow") 
    display(p)
end

function example()
    urdf = "C:\\Users\\jcalbert\\Documents\\GitHub\\Dionysos.jl\\problems\\double_pendulum\\double_pendulum.urdf"
    simulate_double_pendulum(urdf)
    plot_double_pendulum(urdf)
end

# example()

end # module
