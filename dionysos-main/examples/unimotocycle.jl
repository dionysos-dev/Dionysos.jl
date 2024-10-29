# Import necessary modules
include("../src/core/core.jl")

using .CoreControls

# Define the control model for the system
model = @control(
    name="MobileCartMPC",
    problem_type=CustomProblem(),
    system_type=Discrete()  # Discrete-time system
)

# Parameters
@parameter(model, "x_r1", 0.5)  # Reference position 1
@parameter(model, "x_r2", 0.5)
@parameter(model, "N", 20)  # Prediction horizon
@parameter(model, "u1_min", 0.2)
@parameter(model, "u1_max", 2.0)
@parameter(model, "u2_min", -1.0)
@parameter(model, "u2_max", 1.0)

# State Variables (x1, x2 for position, x3 for orientation)
@statevar(model, "x1", Reals())
@statevar(model, "x2", Reals())
@statevar(model, "x3", Reals(), bounds=(-2*pi, 2*pi))  # Modulo operation handled later


# Control Inputs (u1 for linear velocity and u2 for angular velocity)
@inputvar(model, "u1", Reals(), bounds=(u1_min, u1_max))
@inputvar(model, "u2", Reals(), bounds=(u2_min, u2_max))

# Dynamics for each state variable
for t in 0:N-1
    @constraint(model, :(x1[t+1] == x1[t] + u1[t] * cos(x3[t])))
    @constraint(model, :(x2[t+1] == x2[t] + u1[t] * sin(x3[t])))
    @constraint(model, :(x3[t+1] == (x3[t] + u2[t]) % (2 * pi)))  # Modulo operation
end

# State Constraints (Nonlinear)
for t in 0:N-1
    @constraint(model, :(x1[t]^2 - x2[t]^2 <= 4))
    @constraint(model, :(4 * x2[t]^2 - x1[t]^2 <= 16))
end

# Stage Cost
for t in 0:N-1
    @constraint(model, :(stage_cost[t] == 100 * ((x1[t] - x_r1)^2 + (x2[t] - x_r2)^2) + u1[t]^2 + u2[t]^2))
end

# Terminal Cost
@constraint(model, :(terminal_cost == 100 * ((x1[N] - x_r1)^2 + (x2[N] - x_r2)^2)))

# Objective: Minimize stage and terminal costs over the horizon
@minimize(model, :(sum(stage_cost[t] for t in 0:N-1) + terminal_cost))

# Solve the problem
algorithm = @uniformgrid(model) # Default algorithm
solution = solve(model, algorithm; horizon=10)

# Print the solution (Optional)
print(solution)
