__precompile__()

module BipedRobot

using RigidBodyDynamics
using RigidBodyDynamics.Contact
using StaticArrays
using Symbolics
using MechanismGeometries
using LightXML

packagepath() = joinpath(@__DIR__, "..", "deps")
urdfpath() = joinpath(packagepath(), "biped_robot.urdf")

function __init__()
    if !isfile(urdfpath())
        error(
            "Could not find $(urdfpath()). Please run `import Pkg; Pkg.build(\"BipedRobot\")`.",
        )
    end
end

function default_contact_model()
    return SoftContactModel(
        hunt_crossley_hertz(; k = 500e3),
        ViscoelasticCoulombModel(0.8, 20e3, 100.0),
    )
end

function mechanism(;
    symbolic::Bool = false,
    use_urdf::Bool = !symbolic,
    T::Type = symbolic ? Num : Float64,
    add_contact_points = true,
    add_flat_ground = true,
    contactmodel = default_contact_model(),
)

    # Define the robot mechanism
    if use_urdf
        mechanism = RigidBodyDynamics.parse_urdf(
            urdfpath();
            scalar_type = T,
            remove_fixed_tree_joints = false,
        )
        visual = URDFVisuals(urdfpath())
        e = root(visual.xdoc)
        l_l_f = parse(
            Float64,
            attribute(e["link"][end]["visual"][1]["geometry"][1]["sphere"]..., "radius"),
        )

    else
        if symbolic
            inertias =
                @variables m_h m_r_t m_r_l m_r_f m_l_t m_l_l m_l_f I_h I_r_t I_r_l I_r_f I_l_t I_l_l I_l_f positive =
                    true
            lengths =
                @variables l_h_x l_h_y l_h_z l_r_t c_r_t l_r_l c_r_l c_r_f l_l_t c_l_t l_l_f l_r_f c_l_l c_l_f real =
                    true

            gravitational_acceleration = @variables g real = true
        else
            error("Not supported yet, ask Francois")
        end

        params = [inertias..., lengths..., gravitational_acceleration...]
        transpose(params)

        axis = SVector(zero(T), one(T), zero(T)) # axis of rotation for each for the knee and hips joints (around y-axis)
        mechanism = Mechanism(RigidBody{T}("world"); gravity = SVector(zero(T), zero(T), g))
        world = root_body(mechanism) # the fixed 'world' rigid body	

        # Attach hips 
        inertia1 = SpatialInertia(
            CartesianFrame3D("hips_link");
            moment = Symbolics.variables(:I_h, 1:3, 1:3),
            com = SVector(l_h_x, l_h_y, l_h_z),
            mass = m_h,
        )

        body1 = RigidBody(inertia1)
        # Attach the hip link to the world via a floating joint named 'world_to_hips'
        joint1 = Joint("world_to_hips", QuaternionFloating{Num}())
        joint1_to_world = Transform3D(
            frame_before(joint1),
            default_frame(world),
            SVector(zero(T), one(T), one(T)),
        )
        attach!(mechanism, world, body1, joint1; joint_pose = joint1_to_world) # attach the link to the mechanism tree

        # We now repeat the process to the other parts of the robot.

        # Attach the right thigh link to the hips link via a revolute joint named 'r_hips_joint'
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
        attach!(mechanism, body1, body2, joint2; joint_pose = joint2_to_body1)

        # Attach the left thigh link to the hips link via a revolute joint named 'l_hips_joint'
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
        attach!(mechanism, body1, body4, joint4; joint_pose = joint4_to_body1)

        # Attach the right leg link to the right thigh link via a revolute joint named 'r_knee_joint'
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
        attach!(mechanism, body2, body3, joint3; joint_pose = joint3_to_body2)

        # Attach the left leg link to the left thigh link via a revolute joint named 'l_knee_joint'
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
        attach!(mechanism, body4, body5, joint5; joint_pose = joint5_to_body4)

        # Attach the left foot to the left leg link via a revolute joint named 'l_foot_link'
        inertia6 = SpatialInertia(
            CartesianFrame3D("l_foot_link");
            moment = I_l_f * axis * transpose(axis),
            com = SVector(zero(T), zero(T), c_l_f),
            mass = m_l_f,
        )
        body6 = RigidBody(inertia6)
        joint6 = Joint("l_foot_joint", Fixed{Num}())
        joint6_to_body5 = Transform3D(
            frame_before(joint6),
            default_frame(body5),
            SVector(zero(T), zero(T), -l_l_f),
        )
        attach!(mechanism, body5, body6, joint6; joint_pose = joint6_to_body5)

        # Attach the right foot to the right leg link via a revolute joint named 'r_foot_link'
        inertia7 = SpatialInertia(
            CartesianFrame3D("r_foot_link");
            moment = I_r_f * axis * transpose(axis),
            com = SVector(zero(T), zero(T), c_r_f),
            mass = m_r_f,
        )
        body7 = RigidBody(inertia7)
        joint7 = Joint("r_foot_joint", Fixed{Num}())
        joint7_to_body4 = Transform3D(
            frame_before(joint7),
            default_frame(body4),
            SVector(zero(T), zero(T), -l_r_f),
        )
        attach!(mechanism, body4, body7, joint7; joint_pose = joint7_to_body4)
    end

    # Configuration of the contact points
    if !symbolic && add_contact_points && contactmodel != nothing
        for side in (:left, :right)
            foot_link = findbody(mechanism, "$(first(string(side)))_foot_link")
            frame = default_frame(foot_link)
            add_contact_point!(
                foot_link,
                ContactPoint(Point3D(frame, zero(T), zero(T), -l_l_f), contactmodel),
            )
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

function setnominal!(Robotstate::MechanismState)
    return set_configuration!(Robotstate, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
end

include("piracy.jl")

end # module
