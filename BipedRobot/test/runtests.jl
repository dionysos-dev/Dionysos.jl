using BipedRobot
using RigidBodyDynamics
using Test

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

"""
    _include_sandbox(filename)

Include the `filename` in a temporary module that acts as a sandbox. (Ensuring
no constants or functions leak into other files.)

This function was taken from `JuMP/docs/make.jl`.
"""
function _include_sandbox(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

for file in ["mass_matrix_double_pendulum_urdf.jl", "Biped_robot.jl"]
    @testset "$file" begin
        _include_sandbox(file)
    end
end
