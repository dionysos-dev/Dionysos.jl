"""
Visualization tools for robot simulation.

This module depends on MeshCat-related packages and should preferably be loaded
only on the main process, not on lightweight worker processes.
"""
module RSVisualization

using RigidBodyDynamics
using MeshCat
using MeshCatMechanisms
using MechanismGeometries
using StaticArrays

include("RSCore.jl")
using .RSCore

export get_visualization_tool,
    set_visualizer,
    update_visualizer!,
    show_frame!,
    animate_trajectory!,
    set_nominal!,
    set_initialbody!

"""
Construct a simulator and its visualizer.
"""
function get_visualization_tool(;
    robot_urdf = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf"),
)
    rs = RSCore.RobotSimulator(;
        fileName = robot_urdf,
        symbolic = false,
        add_contact_points = true,
        add_gravity = true,
        add_flat_ground = true,
    )
    vis = set_visualizer(; mechanism = rs.mechanism, fileName = robot_urdf)
    return rs, vis
end

"""
Construct the visualizer for a mechanism from a URDF file.
"""
function set_visualizer(; mechanism::Mechanism, fileName::String = "ZMP_2DBipedRobot.urdf")
    urdfpath() = fileName
    return MechanismVisualizer(mechanism, URDFVisuals(urdfpath()))
end

"""
Update the visualizer from the simulator state.
"""
function update_visualizer!(rs::RSCore.RobotSimulator, vis::MechanismVisualizer)
    return set_configuration!(vis, RigidBodyDynamics.configuration(rs.state))
end

"""
Show the default frame of every body in the visualizer.
"""
function show_frame!(rs::RSCore.RobotSimulator, vis::MechanismVisualizer)
    for body in RigidBodyDynamics.bodies(rs.mechanism)
        frame = RigidBodyDynamics.default_frame(body)
        setelement!(vis, frame)
    end
    return nothing
end

"""
Set the simulator to a nominal state and immediately update the visualizer.
"""
function set_nominal!(
    rs::RSCore.RobotSimulator,
    vis::MechanismVisualizer,
    boom::AbstractVector,
    actuators::AbstractVector,
    foot::AbstractVector,
)
    RSCore.set_nominal_state!(rs, boom, actuators, foot)
    return update_visualizer!(rs, vis)
end

"""
Set the simulator to its default initial state and immediately update the visualizer.
"""
function set_initialbody!(rs::RSCore.RobotSimulator, vis::MechanismVisualizer)
    RSCore.set_initialbody_state!(rs)
    return update_visualizer!(rs, vis)
end

"""
Animate a reduced trajectory in the visualizer.
"""
function animate_trajectory!(vis::MechanismVisualizer, x_traj; dt = 0.1)
    fps = max(1, round(Int, 1 / dt))
    anim = MeshCat.Animation(vis.visualizer; fps = fps)

    for (k, x) in enumerate(x_traj)
        q, _ = RSCore.x_to_mechanism_state(x)
        atframe(anim, k) do
            return set_configuration!(vis, q)
        end
    end

    setanimation!(vis, anim)
    return nothing
end

end
