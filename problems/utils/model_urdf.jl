using RigidBodyDynamics, Symbolics
# global urdf_path

function get_mechanism(; wanted_mech = "double_pendulum")
    # Load the URDF file
    if wanted_mech == "double_pendulum"
        urdf_path = "../utils/urdfs/double_pendulum_right/robot.urdf"
    else
        # urdf_path = "../utils/urdfs/doublependulum_new_knee_fixed.urdf"
        urdf_path = "../utils/urdfs/single_pendulum_hip/robot.urdf"
    end
    mechanism = parse_urdf(Float64, urdf_path)
    remove_fixed_tree_joints!(mechanism)

    return mechanism
end

# @variables q1 q2 v1 v2 tau1 tau2 real = true
# q = [q1, q2]
# v = [v1, v2]
# global M, C_G
# M = simplify.(mass_matrix(MechanismState(doublependulum, q, v)));
# C_G = simplify.(dynamics_bias(MechanismState(doublependulum, q, v)));

# MeshCatMechanisms, MeshCat doesn't work if added in the same environment as ModelingToolkit... 
# function animate_single_pendulum(mechanism, qs, ts)
#     mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf_path))
#     animation = Animation(mvis, ts, qs)
#     setanimation!(mvis, animation)

#     # Create a MechanismVisualizer and visualize
#     MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 1.)
# end
