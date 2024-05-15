# using ModelingToolkit, RigidBodyDynamics, DifferentialEquations
# using ModelingToolkit: t_nounits as t, D_nounits as D

include(joinpath("..", "model_urdf.jl"))
include("double_pendulum_system.jl")

mechanism = get_mechanism(;wanted_mech = "double_pendulum")

# Estimate the next state from a given state u0 = [q_1, q̇_1, q_2, q̇_2], for a given input, in a given time tspan
input = [0., .4]
tspan = 5.
stateVectorInit = [0., 0., 0., 0.] # [q_1, q̇_1, q_2, q̇_2]
sys, model = double_pendulum_system.system(mechanism)
prob = double_pendulum_system.problem(sys, constant_input=input)
stateVectorNext = double_pendulum_system.get_next_state(sys, prob, stateVectorInit, input, tspan) # then get the next state

# Once the system is defined, it is no more necessary to do so
input = [1., .2]
tspan = 1.
stateVectorInit = [0., 0., 0., 0.]
stateVectorNext2 = double_pendulum_system.get_next_state(sys, prob, stateVectorInit, input, tspan)

# We can also linearise the system around a given operation point
operation_point = [0., 0., 0., 0.]
matrices = double_pendulum_system.linear_system(sys, model, operation_point)

# matrices hold (; A, B, C, D)