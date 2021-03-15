module DoublePendulum
using RigidBodyDynamics, MeshCat, MeshCatMechanisms, SymPy
using LinearAlgebra,StaticArrays,Random
function test()
    lc1 = -0.5  #center of mass
    l1 = -1.    #length of body 1
    m1 = 1.
    I1 = 0.333 # about joint instead of CoM in URDF
    lc2 = -1.
    l2 = -2.
    m2 = 1.
    I2 = 1.33 # about joint instead of CoM in URDF
    g = -9.81

    axis = SVector(0., 1., 0.)

    world = RigidBody{Float64}("world")
    doublependulum = Mechanism(world; gravity = SVector(0, 0, g))

    # create first body and attach it to the world via a revolute joint
    inertia1 = SpatialInertia(CartesianFrame3D("upper_link"), moment=I1 * axis * axis', com=SVector(0, 0, lc1), mass=m1)
    body1 = RigidBody(inertia1)
    joint1 = Joint("shoulder", Revolute(axis))
    joint1_to_world = one(Transform3D, joint1.frame_before, default_frame(world))
    attach!(doublependulum, world, body1, joint1, joint_pose = joint1_to_world)

    inertia2 = SpatialInertia(CartesianFrame3D("lower_link"), moment=I2 * axis * axis', com=SVector(0, 0, lc2), mass=m2)
    body2 = RigidBody(inertia2)
    joint2 = Joint("elbow", Revolute(axis))
    joint2_to_body1 = Transform3D(joint2.frame_before, default_frame(body1), SVector(0, 0, l1))
    attach!(doublependulum, body1, body2, joint2, joint_pose = joint2_to_body1)

    state = MechanismState(doublependulum)

    # ## The state of a `Mechanism`

    # A `Mechanism` stores the joint/rigid body layout, but no state information. State information is separated out into a `MechanismState` object:

    state = MechanismState(doublependulum)


    # Let's first set the configurations and velocities of the joints:

    set_configuration!(state, joint1, 0.3)
    set_configuration!(state, joint2, 0.4)
    set_velocity!(state, joint1, 1.)
    set_velocity!(state, joint2, 2.)
    # **Important**: a `MechanismState` contains cache variables that depend on the configurations and velocities of the joints. These need to be invalidated when the configurations and velocities are changed. To do this, call
    setdirty!(state)

    ts, qs, vs = simulate(state, 5., Δt = 1e-3)

    q = configuration(state)
    v = velocity(state)


    # a bit simpler double pendulum
    urdf = "doublependulum.urdf"
    vis = MechanismVisualizer(doublependulum, URDFVisuals(urdf))

    t, q, v = simulate(state, 5.0)
    animation = Animation(vis, t, q)
    setanimation!(vis, animation)
    open(vis)

end


function symbolic_model()

    # ## Create symbolic parameters
    # * Masses: $m_1, m_2$
    # * Mass moments of inertia (about center of mass): $I_1, I_2$
    # * Link lengths: $l_1, l_2$
    # * Center of mass locations (w.r.t. preceding joint axis): $c_1, c_2$
    # * Gravitational acceleration: $g$

    inertias = @syms m_1 m_2 I_1 I_2 positive = true
    lengths = @syms l_1 l_2 c_1 c_2 real = true
    gravitational_acceleration = @syms g real = true
    params = [inertias..., lengths..., gravitational_acceleration...]
    transpose(params)


    # ## Create double pendulum `Mechanism`
    # A `Mechanism` contains the joint layout and inertia parameters, but no state information.

    T = Sym # the 'scalar type' of the Mechanism we'll construct
    axis = SVector(zero(T), one(T), zero(T)) # axis of rotation for each of the joints
    double_pendulum = Mechanism(RigidBody{T}("world"); gravity = SVector(zero(T), zero(T), g))
    world = root_body(double_pendulum) # the fixed 'world' rigid body

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

    # Set the joint configuration vector of the MechanismState to a new vector of symbolic variables
    q = configuration(x)
    for i in eachindex(q)
        q[i] = symbols("q_$i", real = true)
    end

    # Set the joint velocity vector of the MechanismState to a new vector of symbolic variables
    v = velocity(x)
    for i in eachindex(v)
        v[i] = symbols("v_$i", real = true)
    end


    # ## Compute dynamical quantities in symbolic form

    # Mass matrix
    M = simplify.(mass_matrix(x))
    println(M)

    # Kinetic energy
    K = simplify(kinetic_energy(x))
    println("Energie cin: ",K)

    # Potential energy
    U = simplify(gravitational_potential_energy(x))
    println("Energie pot: ",U)

    #dynamics!(v,result,state,x)
    #println(result)


end


#control! will be called at each time step of the simulation and allows you to specify joint torques
#given the time and the state of the Mechanism.
#torques is a vector of size equal to the number of joint.
function control!(torques::AbstractVector, t, state::MechanismState)
    #println(torques)
    torques[1] = 0.0
    torques[2] = 0.0
    #rand!(torques) # for example
end

#remark, the acrobot is exactly has the previous pendulum (just the visual is changing)
function acrobot_urdf()
    urdf = "Acrobot.urdf"
    mechanism = parse_urdf(urdf)

    q = [π,Float64(0.1)]
    v = [0.0,0.0]
    state = MechanismState(mechanism, q, v)
    t, q, v = simulate(state, 5.0,control!; Δt = 1e-3)

    q1_min = q[1][1]; q1_max = q[1][1]
    q2_min = q[1][2]; q2_max = q[1][2]
    v1_min = v[1][1]; v1_max = v[1][1]
    v2_min = v[1][2]; v2_max = v[1][2]
    for i=1:length(t)
        q1_min = min(q1_min,q[i][1]); q1_max = max(q1_max,q[i][1])
        q2_min = min(q2_min,q[i][2]); q2_max = max(q2_max,q[i][2])
        v1_min = min(v1_min,v[i][1]); v1_max = max(v1_max,v[i][1])
        v2_min = min(v2_min,v[i][2]); v2_max = max(v2_max,v[i][2])
    end
    println("q1_min : ", q1_min); println("q1_max : ", q1_max)
    println("q2_min : ", q2_min); println("q2_max : ", q2_max)
    println("v1_min : ", v1_min); println("v1_max : ", v1_max)
    println("v2_min : ", v2_min); println("v2_max : ", v2_max)
    #println(q)
    #=#or
    #Let's write a simple controller that just applies $10 \sin(t)$ at the elbow joint and adds some damping at the shoulder joint:
    shoulder, elbow = joints(mechanism)
    function simple_control!(torques::AbstractVector, t, state::MechanismState)
        torques[velocity_range(state, shoulder)] .= -1 .* velocity(state, shoulder)
        torques[velocity_range(state, elbow)] .= 10 * sin(t)
    end
    state = MechanismState(mechanism)
    zero_velocity!(state)
    set_configuration!(state, shoulder, 0.7)
    set_configuration!(state, elbow, -0.8)
    t, q, v = simulate(state, 10.0,simple_control!; Δt = 1e-3)
    =#

    #Simulation
    vis = Visualizer()
    open(vis)  # open the visualizer in a separate tab/window
    mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf),vis)
    animation = Animation(mvis, t, q)
    setanimation!(mvis, animation)

    # example of possibilities for printing
    lower_arm = findbody(mechanism, "lower_link")
    setelement!(mvis, default_frame(lower_arm))
    upper_arm = findbody(mechanism, "upper_link")
    setelement!(mvis, default_frame(upper_arm))
    #radius = 0.05
    #name = "my_point"
    #setelement!(mvis, Point3D(default_frame(lower_arm), 0.2, 0.2, 0.2), radius, name)
    #setelement!(mvis, c, radius, name)
    # upper_arm = findbody(mechanism, "upper_link")
    # point = Point3D(default_frame(upper_arm), 0., 0, -2)
end

println()
#test()
#symbolic_model()
acrobot_urdf()
end










#=
#example de structure pour creer le controller
function make_controller(state::MechanismState,body::RigidBody,point::Point3D)
    mechanism = state.mechanism
    world = root_frame(mechanism)
    joint_path = path(mechanism, root_body(mechanism), body)
    point_world = transform(state, point, root_frame(mechanism))
    Jp = point_jacobian(state, joint_path, point_world)
    v̇ = similar(velocity(state))

    function controller!(τ, t, state)
        desired = Point3D(world, circle_origin .+ radius .* SVector(sin(t / ω), 0, cos(t / ω)))
        point_in_world = transform_to_root(state, body) * point
        point_jacobian!(Jp, state, joint_path, point_in_world)
        Kp = 200
        Kd = 20
        Δp = desired - point_in_world
        v̇ .= Kp * Array(Jp)' * Δp.v .- 20 .* velocity(state)
        τ .= inverse_dynamics(state, v̇)
    end
end
controller! = make_controller(state, body, point, circle_origin, radius, ω)
=#





#=

function four_bar_linkage()
    # ## Model definition

    # We're going to create a [four-bar linkage](https://en.wikipedia.org/wiki/Four-bar_linkage) that looks like this:
    # ![fourbar](fourbar.jpg)
    #
    # We'll 'cut' the mechanism at joint 4: joints 1, 2, and 3 will be part of the spanning tree of the mechanism, but joint 4 will be a 'loop joint' (see e.g. Featherstone's 'Rigid Body Dynamics Algorithms'), for which the dynamics will be enforced using Lagrange multipliers.
    #

    #pattern de presentation (changer le contenu)
    # First, we'll define some relevant constants:

    ## gravitational acceleration
    g = -9.81

    ## link lengths
    l_0 = 1.10
    l_1 = 0.5
    l_2 = 1.20
    l_3 = 0.75

    ## link masses
    m_1 = 0.5
    m_2 = 1.0
    m_3 = 0.75

    ## link center of mass offsets from the preceding joint axes
    c_1 = 0.25
    c_2 = 0.60
    c_3 = 0.375

    ## moments of inertia about the center of mass of each link
    I_1 = 0.333
    I_2 = 0.537
    I_3 = 0.4

    ## Rotation axis: negative y-axis
    axis = SVector(0., -1., 0.);

    # Construct the world rigid body and create a new mechanism:

    world = RigidBody{Float64}("world")
    fourbar = Mechanism(world; gravity = SVector(0., 0., g))

    # Next, we'll construct the spanning tree of the mechanism,
    # consisting of bodies 1, 2, and 3 connected by joints 1, 2, and 3.
    # Note the use of the `moment_about_com` keyword (as opposed to `moment`):

    # Link 1 and joint 1:

    joint1 = Joint("joint1", Revolute(axis))
    inertia1 = SpatialInertia(frame_after(joint1),
        com=SVector(c_1, 0, 0),
        moment_about_com=I_1*axis*transpose(axis),
        mass=m_1)
    link1 = RigidBody(inertia1)
    before_joint1_to_world = one(Transform3D,
        frame_before(joint1), default_frame(world))
    attach!(fourbar, world, link1, joint1,
        joint_pose = before_joint1_to_world)

    # Link 2 and joint 2:

    joint2 = Joint("joint2", Revolute(axis))
    inertia2 = SpatialInertia(frame_after(joint2),
        com=SVector(c_2, 0, 0),
        moment_about_com=I_2*axis*transpose(axis),
        mass=m_2)
    link2 = RigidBody(inertia2)
    before_joint2_to_after_joint1 = Transform3D(
        frame_before(joint2), frame_after(joint1), SVector(l_1, 0., 0.))
    attach!(fourbar, link1, link2, joint2,
        joint_pose = before_joint2_to_after_joint1)

    # Link 3 and joint 3:

    joint3 = Joint("joint3", Revolute(axis))
    inertia3 = SpatialInertia(frame_after(joint3),
        com=SVector(l_0, 0., 0.),
        moment_about_com=I_3*axis*transpose(axis),
        mass=m_3)
    link3 = RigidBody(inertia3)
    before_joint3_to_world = Transform3D(
        frame_before(joint3), default_frame(world), SVector(l_0, 0., 0.))
    attach!(fourbar, world, link3, joint3, joint_pose = before_joint3_to_world)

    # Finally, we'll add joint 4 in almost the same way we did the other joints,
    # with the following exceptions:
    # 1. both `link2` and `link3` are already part of the `Mechanism`, so the `attach!`
    #    function will figure out that `joint4` will be a loop joint.
    # 2. instead of using the default (identity) for the argument that specifies the
    #    transform from the successor of joint 4 (i.e., link 3) to the frame directly after
    # joint 4, we'll specify a transform that incorporates the $l_3$ offset.

    ## joint between link2 and link3
    joint4 = Joint("joint4", Revolute(axis))
    before_joint4_to_joint2 = Transform3D(
        frame_before(joint4), frame_after(joint2), SVector(l_2, 0., 0.))
    joint3_to_after_joint4 = Transform3D(
        frame_after(joint3), frame_after(joint4), SVector(-l_3, 0., 0.))
    attach!(fourbar, link2, link3, joint4,
        joint_pose = before_joint4_to_joint2, successor_pose = joint3_to_after_joint4)

    # Note the additional non-tree joint in the printed `Mechanism` summary.

    # ## Simulation

    # As usual, we'll now construct a `MechanismState` and `DynamicsResult` for the
    # four-bar `Mechanism`. We'll set some initial conditions for a simulation, which
    # were solved for a priori using a nonlinear program (not shown here).

    state = MechanismState(fourbar)
    result = DynamicsResult(fourbar);

    set_configuration!(state, joint1, 1.6707963267948966) # θ
    set_configuration!(state, joint2, -1.4591054166649482) # γ
    set_configuration!(state, joint3, 1.5397303602625536) # ϕ

    set_velocity!(state, joint1, 0.5)
    set_velocity!(state, joint2, -0.47295)
    set_velocity!(state, joint3, 0.341)

    # Next, we'll do a 3-second simulation:

    ts, qs, vs = simulate(state, 3., Δt = 1e-2);

    # ## Visualization

    # For visualization, we'll use [`MeshCatMechanisms`](https://github.com/JuliaRobotics/MeshCatMechanisms.jl),
    # an external package based on RigidBodyDynamics.jl.



    mvis = MechanismVisualizer(fourbar, Skeleton(inertias=false));
    # And animate:
    MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 1.);

end
=#
