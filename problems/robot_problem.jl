module RobotProblem

bipedrobotpackagepath() = joinpath(@__DIR__, "..", "BipedRobot", "src", "BipedRobot.jl")
include(bipedrobotpackagepath())

using RigidBodyDynamics
using RigidBodyDynamics.Contact
using Random
using StaticArrays
using MathematicalSystems
using Dionysos
using .BipedRobot

function vectorFieldBipedRobot(x, u)
    q = x[1:num_positions(state)]
    v = x[(num_positions(state) + 1):(num_positions(state) + num_velocities(state))]
    s = x[(num_positions(state) + num_velocities(state) + 1):end]

    set_configuration!(state, q)
    set_velocity!(state, v)
    set_additional_state!(state, s)

    ẋ = similar(x)

    dynamics!(ẋ, DynamicsResult(state.mechanism), state, x, u)

    return ẋ
end

function system()
    robot = BipedRobot.mechanism(;
        symbolic = false,
        add_contact_points = true,
        add_flat_ground = true,
    )
    state = MechanismState(robot)

    sys = BlackBoxControlContinuousSystem(
        vectorFieldBipedRobot,
        num_positions(state) + num_velocities(state) + num_additional_states(state),
        num_velocities(state),
    )

    return sys
end

# TODO: create the Dionysos.Problem.OptimalControlProblem

end
