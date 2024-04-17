# Example of use

# We first include the urdf and the system definition
include("single_pendulum_urdf.jl")
include("single_pendulum_system.jl")

# Define the single-pendulum mechanism (RigidBodyDynamics)
mechanism = single_pendulum_urdf.single_pendulum_mechanism()

# Estimate the next state from a given state u0 = [q, q̇], for a given input, in a given time tspan
input = 2.
tspan = 5.
u0 = [0., 0.]
sys, model = single_pendulum_system.system(mechanism) # First define the system
next_u = single_pendulum_system.get_next_state(sys, u0, input, tspan) # then get the next state

# Once the system is defined, it is no more necessary to do so
input = 1.
tspan = 1.
u0 = [0., 0.]
next_u2 = single_pendulum_system.get_next_state(sys, u0, input, tspan)


# We can also linearise the system around a given operation point
operation_point = [0., 0.]
matrices = single_pendulum_system.linear_system(sys, model, operation_point)

# matrices hold (; A, B, C, D)
