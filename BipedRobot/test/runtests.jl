using BipedRobot
using RigidBodyDynamics
using Test

# Test for mechanism via URDF + setnominal :
#mechanism_urdf = BipedRobot.mechanism(symbolic = false, add_contact_points = true, add_flat_ground = true)
#state_urdf = MechanismState(mechanism_urdf)
#BipedRobot.setnominal!(state_urdf)

# Otherwise ambiguous normalize
normalize(v::FreeVector3D, p::Real) = RigidBodyDynamics.Spatial.normalize(v, p) 

# Test for mechanism via symbolic variables :
mechanism_sym = BipedRobot.mechanism(;
    symbolic = false,
    add_contact_points = true,
    add_flat_ground = true,
)
state_sym = MechanismState(mechanism_sym)
M = mass_matrix(state_sym)
@test size(M) == (10, 10)

mechanism_sym = BipedRobot.mechanism(;
    symbolic = true,
    add_contact_points = true,
    add_flat_ground = true,
)
state_sym = MechanismState(mechanism_sym)
mass_matrix(state_sym)
M = BipedRobot.simplify.(mass_matrix(state_sym))
@test size(M) == (10, 10)
