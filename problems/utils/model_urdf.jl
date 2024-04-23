using RigidBodyDynamics

function get_mechanism(; wanted_mech = "double_pendulum")
    # Load the URDF file
    if wanted_mech == "double_pendulum"
        urdf_path = joinpath(@__DIR__, "urdfs", "double_pendulum_right", "robot.urdf")
    else
        urdf_path = joinpath(@__DIR__, "urdfs", "single_pendulum_hip", "robot.urdf")
    end
    mechanism = parse_urdf(Float64, urdf_path)
    remove_fixed_tree_joints!(mechanism)

    return mechanism
end
