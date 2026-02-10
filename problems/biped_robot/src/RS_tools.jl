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
    fileName::String = "ZMP_2DBipedRobot.urdf",
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
    return RobotSimulator(robot, state)#, Array[], Array[], Array[])
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
    fileName::String = "ZMP_bipedalRobot.urdf",
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

end
