using BipedRobot
using RigidBodyDynamics

# Test for mechanism via URDF + setnominal :
#mechanism_urdf = BipedRobot.mechanism(symbolic = false, add_contact_points = true, add_flat_ground = true)
#state_urdf = MechanismState(mechanism_urdf)
#BipedRobot.setnominal!(state_urdf)

# Test for mechanism via symbolic variables :
mechanism_sym = BipedRobot.mechanism(;
    symbolic = false,
    add_contact_points = true,
    add_flat_ground = true,
)
state_sym = MechanismState(mechanism_sym)
simplify.(mass_matrix(state_sym))
