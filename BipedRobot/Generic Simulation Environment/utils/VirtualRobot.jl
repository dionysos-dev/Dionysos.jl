mutable struct VirtualRobot
    mechanism::Mechanism
    state::MechanismState
    vis::MechanismVisualizer
    contacts_points::Vector{Array}
    CoM::Array
    torques::Array
end # End Module 

function VirtualRobot(;
    fileName::String = "ZMP_2DBipedRobot.urdf",
    contactmodel = default_contact_model(),
    contacts_points::Vector{Array} = [],
    add_contact_points = true,
    add_flat_ground = true,
    add_gravity = true,
)
    urdfpath() = joinpath(packagepath(), fileName)
    println(urdfpath())
    mechanism = RigidBodyDynamics.parse_urdf(urdfpath())
    visual = MechanismVisualizer(mechanism, URDFVisuals(urdfpath()))

    if add_gravity
        mechanism.gravitational_acceleration =
            FreeVector3D(root_frame(mechanism), 0, 0, -9.81)
        println("Gravity added")
    else
        mechanism.gravitational_acceleration = FreeVector3D(root_frame(mechanism), 0, 0, 0)
    end
    if add_contact_points && contactmodel !== nothing
        for side in (:left, :right)
            foot_link = findbody(mechanism, "$(first(string(side)))_foot_link")
            frame = default_frame(foot_link)
            for (i, cp) in enumerate(contacts_points)
                point = Point3D(frame, cp[1], cp[2], cp[3])
                add_contact_point!(foot_link, ContactPoint(point, contactmodel))
                # Show contacts points
                setelement!(visual, point, 0.005, "cp$(i)")
            end
        end
        println("Contatcs added")
    end
    if add_flat_ground
        frame = root_frame(mechanism)
        ground =
            HalfSpace3D(Point3D(frame, 0.0, 0.0, 0.0), FreeVector3D(frame, 0.0, 0.0, 1.0))
        add_environment_primitive!(mechanism, ground)
        println("Ground added")
    end
    state = MechanismState(mechanism)
    return VirtualRobot(mechanism, state, visual, contacts_points, [], [])
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

Set the actual robot `mechansim` to his nominal state 
"""
function set_nominal!(vr::VirtualRobot, config::Vector)
    set_configuration!(vr.state, config)
    zero_velocity!(vr.state)
    return update_visulizer!(vr)
end

"""

Set the actual robot `mechansim` to his initial state 
"""
function set_initialbody!(vr::VirtualRobot)
    config = zeros(length(vr.state.v))
    set_configuration!(vr.state, config)
    zero_velocity!(vr.state)
    return update_visulizer!(vr)
end

"""

Show the joint frame in the visualiser 
"""
function show_frame!(vr::VirtualRobot)
    # Show the frame of each bodies
    robot_bodies = RigidBodyDynamics.bodies(vr.mechanism)
    for body in robot_bodies
        frame = RigidBodyDynamics.default_frame(body)
        setelement!(vr.vis, frame)
    end
end

"""

Update the visulatiser to the new state 
"""
function update_visulizer!(vr::VirtualRobot)
    return set_configuration!(vr.vis, RigidBodyDynamics.configuration(vr.state))
end
