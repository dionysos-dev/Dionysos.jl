"""
A Robot Simulator allows to define the simulation envrionnment.

A `RobotSimulator` has 5 differents outputs :
    1.  `mechanism` which is a data structure which contains the state of the robot and more 
    2.  `state` which is the actual state of the robot 
"""
module RS_tools
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using StructArrays
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using MeshCat, MeshCatMechanisms, Blink
using MechanismGeometries
using LightXML
using GeometryTypes
using StaticArrays

export RobotSimulator,
    getMechanism, set_nominal!, set_initialbody!, update_visulizer!, show_frame!, simulate

struct RobotSimulator
    mechanism::Mechanism
    state::MechanismState
end

"""
Constructor 
"""
function RobotSimulator(;
    fileName::String = "",
    symbolic::Bool = false,
    use_urdf::Bool = !symbolic,
    T::Type = symbolic ? Num : Float64,
    add_contact_points = true,
    add_flat_ground = true,
    add_gravity = true,
    contactmodel = default_contact_model(),
)
    robot = getMechanism(;
        fileName,
        symbolic,
        use_urdf,
        T,
        add_contact_points,
        add_flat_ground,
        add_gravity,
        contactmodel,
    )
    state = MechanismState(robot)
    return RobotSimulator(robot, state)
end

"""
Define the contact model and their parameters 
"""
function default_contact_model()
    return SoftContactModel(
        hunt_crossley_hertz(; k = 500e2),
        ViscoelasticCoulombModel(0.8, 20e3, 100.0),
    )
end

"""
Return the URDF  `mechanism` of the robot. 

Describe the simulation environement as well. 
"""
function getMechanism(;
    fileName::String = "",
    symbolic::Bool = false,
    use_urdf::Bool = !symbolic,
    T::Type = symbolic ? Num : Float64,
    add_contact_points = true,
    add_flat_ground = true,
    add_gravity = true,
    contactmodel = default_contact_model(),
)
    urdfpath() = fileName

    # Define the robot mechanism
    if use_urdf
        mechanism = RigidBodyDynamics.parse_urdf(urdfpath())
        visual = URDFVisuals(urdfpath())
        if (add_gravity)
            mechanism.gravitational_acceleration =
                FreeVector3D(root_frame(mechanism), 0, 0, -9.81)
        else
            mechanism.gravitational_acceleration =
                FreeVector3D(root_frame(mechanism), 0, 0, 0)
        end
    else
        error("Not implemented yet ! Please use URDF file !")
    end
    # Configuration of the contact points
    if !symbolic && add_contact_points && contactmodel !== nothing
        e = root(visual.xdoc)
        s = attribute(e["link"][end]["visual"][1]["geometry"][1]["box"]..., "size")

        numbers = split(s, " ")
        left_foot_length = parse(Float64, numbers[1])
        left_foot_width = parse(Float64, numbers[2])
        left_foot_height = parse(Float64, numbers[3])

        for side in (:left, :right)
            foot_link = findbody(mechanism, "$(first(string(side)))_foot_link")
            frame = default_frame(foot_link)
            cp0 = Point3D(frame, zero(T), zero(T), -left_foot_height)
            add_contact_point!(foot_link, ContactPoint(cp0, contactmodel))
            for div in 2:2
                for sign in [-1, 1]
                    cp = Point3D(
                        frame,
                        sign * left_foot_length / div,
                        sign * left_foot_width / div,
                        -left_foot_height,
                    )
                    add_contact_point!(foot_link, ContactPoint(cp, contactmodel))
                    cp = Point3D(
                        frame,
                        -sign * left_foot_length / div,
                        sign * left_foot_width / div,
                        -left_foot_height,
                    )
                    add_contact_point!(foot_link, ContactPoint(cp, contactmodel))
                end
            end
        end
    end
    # Configuration of the ground
    if !symbolic && add_flat_ground
        frame = root_frame(mechanism)
        ground =
            HalfSpace3D(Point3D(frame, 0.0, 0.0, 0.0), FreeVector3D(frame, 0.0, 0.0, 1.0))
        add_environment_primitive!(mechanism, ground)
    end
    return mechanism
end

function get_visualization_tool(;
    robot_urdf = joinpath(@__DIR__, "..", "deps/ZMP_2DBipedRobot_nodamping.urdf"),
)
    rs = RS_tools.RobotSimulator(;
        fileName = robot_urdf,
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    )
    vis = RS_tools.set_visulalizer(; mechanism = rs.mechanism, fileName = robot_urdf)
    return rs, vis
end

"""
Set the actual robot `mechansim` to his nominal state 
"""
function set_nominal!(
    rs::RobotSimulator,
    vis::MechanismVisualizer,
    boom::VecOrMat,
    actuators::VecOrMat,
    foot::VecOrMat,
)
    config = vec(reduce(hcat, [boom; actuators; foot]))
    set_configuration!(rs.state, config)
    zero_velocity!(rs.state)
    set_configuration!(vis, configuration(rs.state))
    return update_visulizer!(rs, vis)
end

"""
Set the actual robot `mechansim` to his initial state 
"""
function set_initialbody!(rs::RobotSimulator, vis::MechanismVisualizer)
    boom = [0, 0]
    foot = [0, 0]
    actuators = [0, 0, 0, 0]
    config = vec(reduce(hcat, [boom; actuators; foot]))
    set_configuration!(rs.state, config)
    zero_velocity!(rs.state)
    return update_visulizer!(rs, vis)
end

"""
Return the visualiser 
"""
function set_visulalizer(; mechanism::Mechanism, fileName::String = "ZMP_2DBipedRobot.urdf")
    urdfpath() = fileName
    vis = MechanismVisualizer(mechanism, URDFVisuals(urdfpath()))
    return vis
end

"""
Update the visualiser to the new state 
"""
function update_visulizer!(rs::RobotSimulator, vis::MechanismVisualizer)
    return set_configuration!(vis, RigidBodyDynamics.configuration(rs.state))
end

"""
Show the joint frame in the visualiser 
"""
function show_frame!(rs::RobotSimulator, vis::MechanismVisualizer)
    # Show the frame of each bodies
    robot_bodies = RigidBodyDynamics.bodies(rs.mechanism)
    for body in robot_bodies
        frame = RigidBodyDynamics.default_frame(body)
        setelement!(vis, frame)
    end
end

function x_to_mechanism_state(x::SVector{6, <:Real})
    # Same constants as in `system()`
    Lthigh = 0.20125
    Lleg = 0.172
    Init_offset = -0.0006559432

    # Build the “partial” joint vectors (8 DoFs here)
    q = SVector{8, Float64}(0.0, 0.0, x[1], x[2], x[3], 0.0, 0.0, 0.0)
    qd = SVector{8, Float64}(0.0, 0.0, x[4], x[5], x[6], 0.0, 0.0, 0.0)

    # Heights (double pendulum)
    zl = Lthigh * cos(q[3]) + Lleg * cos(q[5] + q[3])
    zr = Lthigh * cos(q[4]) + Lleg * cos(q[6] + q[4])

    # Boom z
    q2 = max(zl, zr) - Lthigh - Lleg + Init_offset
    q = SVector{8, Float64}(q[1], q2, q[3], q[4], q[5], q[6], q[7], q[8])

    # Identify contact leg for boom x velocity terms
    i1 = zl > zr ? 3 : 4
    i2 = zl > zr ? 5 : 6

    xboom = Lthigh * sin(q[i1]) + Lleg * sin(q[i2] + q[i1])
    ẋboom = Lthigh * qd[i1] * cos(q[i1]) + Lleg * (qd[i1] + qd[i2]) * cos(q[i1] + q[i2])
    żboom = -(Lthigh * qd[i1] * sin(q[i1]) + Lleg * (qd[i1] + qd[i2]) * sin(q[i1] + q[i2]))

    q1 = xboom
    qd1 = ẋboom
    qd2 = żboom

    # Keep feet mostly horizontal
    qd7 = -(qd[3] + qd[5])
    qd8 = -(qd[4] + qd[6])

    # Final full configuration and velocity (8D)
    q = SVector{8, Float64}(q1, q2, q[3], q[4], q[5], q[6], -(q[3] + q[5]), -(q[4] + q[6]))
    qd = SVector{8, Float64}(qd1, qd2, qd[3], qd[4], qd[5], qd[6], qd7, qd8)

    return q, qd
end

function animate_trajectory!(vis::MechanismVisualizer, x_traj; dt = 0.1)
    fps = max(1, round(Int, 1 / dt))
    anim = MeshCat.Animation(vis.visualizer; fps = fps)

    for (k, x) in enumerate(x_traj)
        q, _ = x_to_mechanism_state(x)
        atframe(anim, k) do
            return set_configuration!(vis, q)
        end
    end

    setanimation!(vis, anim)
    return nothing
end

end
