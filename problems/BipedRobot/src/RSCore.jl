"""
Core robot-simulation tools without visualization dependencies.

This module is intended for computation-heavy or distributed workloads.
It avoids loading MeshCat / Blink / visualization packages on workers.
"""
module RSCore

using RigidBodyDynamics
using RigidBodyDynamics.Contact
using LightXML
using StaticArrays

export RobotSimulator,
    default_contact_model,
    getMechanism,
    get_visual_geometry_dimensions,
    set_nominal_state!,
    set_initialbody_state!,
    x_to_mechanism_state

"""
A lightweight robot simulator container.

Fields:
- `mechanism`: robot mechanism and environment
- `state`: current mechanism state
"""
struct RobotSimulator
    mechanism::Mechanism
    state::MechanismState
end

"""
Default soft contact model.
"""
function default_contact_model()
    return SoftContactModel(
        hunt_crossley_hertz(; k = 500e2),
        ViscoelasticCoulombModel(0.8, 20e3, 100.0),
    )
end

"""
Read the first box geometry dimensions from a URDF visual tag.

Returns `(length, width, height)`.
"""
function get_visual_geometry_dimensions(fileName::String)
    xdoc = parse_file(fileName)
    rootnode = root(xdoc)

    link_nodes = get_elements_by_tagname(rootnode, "link")
    isempty(link_nodes) && error("No <link> found in URDF file: $fileName")

    last_link = link_nodes[end]
    visual_nodes = get_elements_by_tagname(last_link, "visual")
    isempty(visual_nodes) && error("No <visual> found in last link of URDF file: $fileName")

    geometry_nodes = get_elements_by_tagname(visual_nodes[1], "geometry")
    isempty(geometry_nodes) &&
        error("No <geometry> found in URDF visual of file: $fileName")

    box_nodes = get_elements_by_tagname(geometry_nodes[1], "box")
    isempty(box_nodes) && error("No <box> found in URDF visual geometry of file: $fileName")

    size_attr = attribute(box_nodes[1], "size")
    size_attr === nothing &&
        error("No size attribute found in URDF box geometry of file: $fileName")

    numbers = split(size_attr, " ")
    length(numbers) == 3 || error("Expected 3 box dimensions in URDF, got: $size_attr")

    return parse(Float64, numbers[1]),
    parse(Float64, numbers[2]),
    parse(Float64, numbers[3])
end

"""
Construct the robot mechanism and optionally augment it with contact points,
flat ground, and gravity.
"""
function getMechanism(;
    fileName::String = "",
    symbolic::Bool = false,
    use_urdf::Bool = !symbolic,
    T::Type = symbolic ? Num : Float64,
    add_contact_points::Bool = true,
    add_flat_ground::Bool = true,
    add_gravity::Bool = true,
    contactmodel = default_contact_model(),
)
    urdfpath() = fileName

    if use_urdf
        mechanism = RigidBodyDynamics.parse_urdf(urdfpath())

        if add_gravity
            mechanism.gravitational_acceleration =
                FreeVector3D(root_frame(mechanism), 0, 0, -9.81)
        else
            mechanism.gravitational_acceleration =
                FreeVector3D(root_frame(mechanism), 0, 0, 0)
        end
    else
        error("Not implemented yet. Please use a URDF file.")
    end

    if !symbolic && add_contact_points && contactmodel !== nothing
        left_foot_length, left_foot_width, left_foot_height =
            get_visual_geometry_dimensions(urdfpath())

        for side in (:left, :right)
            foot_link = findbody(mechanism, "$(first(string(side)))_foot_link")
            frame = default_frame(foot_link)

            cp0 = Point3D(frame, zero(T), zero(T), -left_foot_height)
            add_contact_point!(foot_link, ContactPoint(cp0, contactmodel))

            for div in 2:2
                for sign in (-1, 1)
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

    if !symbolic && add_flat_ground
        frame = root_frame(mechanism)
        ground =
            HalfSpace3D(Point3D(frame, 0.0, 0.0, 0.0), FreeVector3D(frame, 0.0, 0.0, 1.0))
        add_environment_primitive!(mechanism, ground)
    end

    return mechanism
end

"""
Constructor for `RobotSimulator`.
"""
function RobotSimulator(;
    fileName::String = "",
    symbolic::Bool = false,
    use_urdf::Bool = !symbolic,
    T::Type = symbolic ? Num : Float64,
    add_contact_points::Bool = true,
    add_flat_ground::Bool = true,
    add_gravity::Bool = true,
    contactmodel = default_contact_model(),
)
    mechanism = getMechanism(;
        fileName,
        symbolic,
        use_urdf,
        T,
        add_contact_points,
        add_flat_ground,
        add_gravity,
        contactmodel,
    )
    state = MechanismState(mechanism)
    return RobotSimulator(mechanism, state)
end

"""
Set the robot state to a nominal configuration.
"""
function set_nominal_state!(
    rs::RobotSimulator,
    boom::AbstractVector,
    actuators::AbstractVector,
    foot::AbstractVector,
)
    config = vec(reduce(hcat, [boom; actuators; foot]))
    set_configuration!(rs.state, config)
    zero_velocity!(rs.state)
    return rs
end

"""
Set the robot state to the default initial body configuration.
"""
function set_initialbody_state!(rs::RobotSimulator)
    boom = [0, 0]
    foot = [0, 0]
    actuators = [0, 0, 0, 0]
    config = vec(reduce(hcat, [boom; actuators; foot]))
    set_configuration!(rs.state, config)
    zero_velocity!(rs.state)
    return rs
end

"""
Convert a reduced 6D state into full mechanism configuration and velocity.
"""
function x_to_mechanism_state(x::SVector{6, <:Real})
    Lthigh = 0.20125
    Lleg = 0.172
    Init_offset = -0.0006559432

    q = SVector{8, Float64}(0.0, 0.0, x[1], x[2], x[3], 0.0, 0.0, 0.0)
    qd = SVector{8, Float64}(0.0, 0.0, x[4], x[5], x[6], 0.0, 0.0, 0.0)

    zl = Lthigh * cos(q[3]) + Lleg * cos(q[5] + q[3])
    zr = Lthigh * cos(q[4]) + Lleg * cos(q[6] + q[4])

    q2 = max(zl, zr) - Lthigh - Lleg + Init_offset
    q = SVector{8, Float64}(q[1], q2, q[3], q[4], q[5], q[6], q[7], q[8])

    i1 = zl > zr ? 3 : 4
    i2 = zl > zr ? 5 : 6

    xboom = Lthigh * sin(q[i1]) + Lleg * sin(q[i2] + q[i1])
    ẋboom = Lthigh * qd[i1] * cos(q[i1]) + Lleg * (qd[i1] + qd[i2]) * cos(q[i1] + q[i2])
    żboom = -(Lthigh * qd[i1] * sin(q[i1]) + Lleg * (qd[i1] + qd[i2]) * sin(q[i1] + q[i2]))

    q1 = xboom
    qd1 = ẋboom
    qd2 = żboom

    qd7 = -(qd[3] + qd[5])
    qd8 = -(qd[4] + qd[6])

    q = SVector{8, Float64}(q1, q2, q[3], q[4], q[5], q[6], -(q[3] + q[5]), -(q[4] + q[6]))
    qd = SVector{8, Float64}(qd1, qd2, qd[3], qd[4], qd[5], qd[6], qd7, qd8)

    return q, qd
end

end
