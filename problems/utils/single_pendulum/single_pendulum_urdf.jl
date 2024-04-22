module single_pendulum_urdf

using RigidBodyDynamics#, Symbolics

function single_pendulum_mechanism()
    # Load the URDF file
    # global urdf_path = "../utils/urdfs/doublependulum_new_knee_fixed.urdf"
    urdf_path = joinpath("..", "utils", "urdfs", "single_pendulum_hip", "robot.urdf")
    mechanism = parse_urdf(Float64, urdf_path)
    remove_fixed_tree_joints!(mechanism)

    return mechanism
end

# MeshCatMechanisms, MeshCat doesn't work if added in the same environment as ModelingToolkit.
# function animate_single_pendulum(mechanism, qs, ts)
#     mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf_path))
#     animation = Animation(mvis, ts, qs)
#     setanimation!(mvis, animation)

#     # Create a MechanismVisualizer and visualize
#     MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 1.)
# end

end