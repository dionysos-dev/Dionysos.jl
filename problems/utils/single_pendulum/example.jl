# Example of use

# We first include the urdf and the system definition
include(joinpath("..", "model_urdf.jl"))
include("single_pendulum_system.jl")

# Define the single-pendulum mechanism (RigidBodyDynamics)
mechanism = get_mechanism(;wanted_mech = "single_pendulum")

# Estimate the next state from a given state u0 = [q, q̇], for a given input, in a given time tspan
input = 2.
tspan = 5.
stateVectorInit = [0., 0.]
# First define the system
sys, model = single_pendulum_system.system(mechanism)
# Then define the problem
prob = single_pendulum_system.problem(sys; u0 = stateVectorInit, tspan = tspan, constant_input = input)
# Then get the next state for a given input and state
stateVectorNext = single_pendulum_system.get_next_state(sys, prob, stateVectorInit, input, tspan)

# Once the problem is defined, it is no more necessary to do so
input = 1.
tspan = 1.
stateVectorInit = [0., 0.]
stateVectorNext2 = single_pendulum_system.get_next_state(sys, prob, stateVectorInit, input, tspan)

# We can also linearise the system around a given operation point
operation_point = [0., 0.]
matrices = single_pendulum_system.linear_system(sys, model, operation_point)

# matrices hold (; A, B, C, D)
