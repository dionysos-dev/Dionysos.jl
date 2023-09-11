using BipedRobot
using RigidBodyDynamics
using Test

# Test for mechanism via URDF + setnominal :
#mechanism_urdf = BipedRobot.mechanism(symbolic = false, add_contact_points = true, add_flat_ground = true)
#state_urdf = MechanismState(mechanism_urdf)
#BipedRobot.setnominal!(state_urdf)

# To avoid ambiguities
using LinearAlgebra, Symbolics, Quaternions
LinearAlgebra.normalize(v::FreeVector3D, p::Real) = FreeVector3D(v.frame, normalize(v.v, p))
Base.:/(q::Quaternion, x::Symbolics.Num) =
    Quaternion(q.s / x.val, q.v1 / x.val, q.v2 / x.val, q.v3 / x.val)

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
