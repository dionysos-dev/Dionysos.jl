# main.jl

# Import necessary modules
include("../src/core/core.jl")

using .CoreControls

using Ipopt
using Plots

# Create a Control model with problem and system types
model = @control(
    name="AutonomousVehiclePathPlanning",
    problem_type=Reachability(),
    system_type=Auto()
)

# Parameters
@parameter(model, "L", 2.5)
@parameter(model, "v_max", 30.0)
@parameter(model, "delta_max", π / 4)
@parameter(model, "Δt", 0.1)
@parameter(model, "N", 10)

# State Variables
@statevar(model, "x", Reals(), bounds=(0, 100))
@statevar(model, "y", Reals(), bounds=(0, 100))
@statevar(model, "psi", Reals(), bounds=(-π, π))

# Control Inputs
@inputvar(model, "v", Reals(), bounds=(0, v_max))
@inputvar(model, "delta", Reals(), bounds=(-delta_max, delta_max))

# Dynamics using dot() for continuous-time system
@constraint(model, :(dot(x, t, Δt, model) == v * cos(psi)))
@constraint(model, :(dot(y, t, Δt, model) == v * sin(psi)))
@constraint(model, :(dot(psi, t, Δt, model) == (v / L) * tan(delta)))

# Initial Conditions
x_start = 0.0
y_start = 0.0
psi_start = 0.0
@constraint(model, :(x_0 == x_start))
@constraint(model, :(y_0 == y_start))
@constraint(model, :(psi_0 == psi_start))

# Terminal Conditions
x_goal = 100.0
y_goal = 100.0
psi_goal = 0.0
@constraint(model, :(x_$N == $x_goal))
@constraint(model, :(y_$N == $y_goal))
@constraint(model, :(psi_$N == $psi_goal))

# Obstacle (Circular)
function add_circular_obstacle(model::Control, x::Variable, y::Variable, x_c::Float64, y_c::Float64, r::Float64, N::Int)
    for t in 0:N
        @constraint(model, (x[t+1] - x_c)^2 + (y[t+1] - y_c)^2 >= r^2)
    end
end

# Obstacle parameters
x_c = 50.0
y_c = 50.0
r = 10.0
add_circular_obstacle(model, x, y, x_c, y_c, r, N)

# Objective Function (Minimize total time)
@minimize(model, :(sum(1 for t in 0:$(N-1))))

# Print the model to verify
print(model)

# Define the algorithm
algorithm = @uniformgrid(model)

# Solve the problem using the specified algorithm
solution = solve(model, algorithm; horizon=N, Δt=Δt)

# Visualize the path
#@visualize(solution, plot="path", labels=["x", "y"])
