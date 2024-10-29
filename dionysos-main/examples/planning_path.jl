# Import necessary modules
include("../src/core/core.jl")

using .CoreControls

# Create a Control model for the vehicle
model = @control(
    name="VehicleMazeNavigation",
    problem_type=Reachability(),
    system_type=Continuous()  # Continuous system type
)

# Parameters
@parameter(model, "u_min", -1.0)
@parameter(model, "u_max", 1.0)
@parameter(model, "alpha", 0.0)  # Initialize alpha, computed later

# State Variables (x1, x2 for position and x3 for orientation)
@statevar(model, "x1", Reals())
@statevar(model, "x2", Reals())
@statevar(model, "x3", Reals())

# Control Inputs (u1 for velocity and u2 for steering angle)
@inputvar(model, "u1", Reals(), bounds=(u_min, u_max))
@inputvar(model, "u2", Reals(), bounds=(u_min, u_max))

# Dynamics for each state variable
@constraint(model, :(dot(x1, t, Δt, model) == u1 * cos(alpha + x3) * cos(alpha)^-1))
@constraint(model, :(dot(x2, t, Δt, model) == u1 * sin(alpha + x3) * cos(alpha)^-1))
@constraint(model, :(dot(x3, t, Δt, model) == u1 * tan(u2)))

# Alpha is dependent on u2
@constraint(model, :(alpha == atan(tan(u2) / 2)))

# Initial Conditions (Specify vehicle's starting position and orientation)
x1_start = 0.0
x2_start = 0.0
x3_start = 0.0
@constraint(model, :(x1_0 == x1_start))
@constraint(model, :(x2_0 == x2_start))
@constraint(model, :(x3_0 == x3_start))

# Terminal Conditions (Target position)
x1_goal = 10.0  # Example target position
x2_goal = 10.0
x3_goal = 0.0
@constraint(model, :(x1_$N == x1_goal))
@constraint(model, :(x2_$N == x2_goal))
@constraint(model, :(x3_$N == x3_goal))

# Obstacle boundaries (provided)
x1_lb = [1.0, 2.2, 2.2, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3]
x1_ub = [1.2, 2.4, 2.4, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0]
x2_lb = [0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2]
x2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6]

# Function to add rectangular obstacle avoidance constraints
function add_rectangular_obstacles(model::Control, x1::Variable, x2::Variable, x1_lb::Vector{Float64}, x1_ub::Vector{Float64}, x2_lb::Vector{Float64}, x2_ub::Vector{Float64}, N::Int)
    num_obstacles = length(x1_lb)
    for i in 1:num_obstacles
        for t in 0:N-1
            # Add constraints to ensure vehicle is outside the obstacle at time t
            @constraint(model, :(x1[t+1] <= $x1_lb[$i] || x1[t+1] >= $x1_ub[$i]))
            @constraint(model, :(x2[t+1] <= $x2_lb[$i] || x2[t+1] >= $x2_ub[$i]))
        end
    end
end

# Add rectangular obstacles to the model
add_rectangular_obstacles(model, model.variables["x1"], model.variables["x2"], x1_lb, x1_ub, x2_lb, x2_ub, N)

# Objective: Minimize time to reach the target position
@minimize(model, :(sum(1 for t in 0:$(N-1))))

# Define the algorithm (e.g., uniform grid)
algorithm = @uniformgrid(model)

# Solve the problem
solution = solve(model, algorithm; horizon=10, Δt=0.1)

# Print the solution (Optional)
##print(solution)

# Visualize the path (Optional)
##Visualize(solution, plot="path", labels=["x1", "x2"])
