# Import necessary modules
include(joinpath(@__DIR__, "../src/core/core.jl"))

#using .CoreControls

# Create a Control model for the vehicle
model = CoreControls.@control(
    name="VehicleMazeNavigation",
    problem_type=CoreControls.Reachability(),
    system_type=CoreControls.Continuous(),  # Continuous system type
)

# Parameters
CoreControls.@parameter(model, "u_min", -1.0)
CoreControls.@parameter(model, "u_max", 1.0)
CoreControls.@parameter(model, "alpha", 0.0)  # Initialize alpha, computed later
CoreControls.@parameter(model, "N", 10)  # Time step

# State Variables (x1, x2 for position and x3 for orientation)
CoreControls.@statevar(model, "x1", CoreControls.Reals())
CoreControls.@statevar(model, "x2", CoreControls.Reals())
CoreControls.@statevar(model, "x3", CoreControls.Reals())

# Control Inputs (u1 for velocity and u2 for steering angle)
CoreControls.@inputvar(model, "u1", CoreControls.Reals(), bounds=(u_min, u_max))
CoreControls.@inputvar(model, "u2", CoreControls.Reals(), bounds=(u_min, u_max))

# Alpha is dependent on u2
CoreControls.@constraint(model, dot(x1, t, Δt, model) == u1 * cos(atan(tan(u2) / 2) + x3) * cos(atan(tan(u2) / 2))^-1)
CoreControls.@constraint(model, dot(x2, t, Δt, model) == u1 * sin(atan(tan(u2) / 2) + x3) * cos(atan(tan(u2) / 2))^-1)
CoreControls.@constraint(model, dot(x3, t, Δt, model) == u1 * tan(u2))

# Initial Conditions (Specify vehicle's starting position and orientation)
CoreControls.@parameter(model, "x1_start", 0.0)
CoreControls.@parameter(model, "x2_start", 0.0)
CoreControls.@parameter(model, "x3_start", 0.0)
CoreControls.@constraint(model, x1_0 == x1_start)
CoreControls.@constraint(model, x2_0 == x2_start)
CoreControls.@constraint(model, x3_0 == x3_start)

# Terminal Conditions (Target position)
CoreControls.@constraint(model, x1_N == 10.0)
CoreControls.@constraint(model, x2_N == 10.0)
CoreControls.@constraint(model, x3_N == 0.0)

# Obstacle boundaries (provided)
x1_lb = [1.0, 2.2, 2.2, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3]
x1_ub = [1.2, 2.4, 2.4, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0]
x2_lb = [0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2]
x2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6]

# Function to add rectangular obstacle avoidance constraints
function add_rectangular_obstacles(model::CoreControls.Control, x1::CoreControls.Variable, x2::CoreControls.Variable, x1_lb::Vector{Float64}, x1_ub::Vector{Float64}, x2_lb::Vector{Float64}, x2_ub::Vector{Float64}, N::Int)
    num_obstacles = length(x1_lb)
    for i in 1:num_obstacles
        CoreControls.@parameter(model, "u_max", 1.0)
        CoreControls.@constraint(
            model,
            x1 <= $(x1_lb[i]) || x1 >= $(x1_ub[i]) ||
            x2 <= $(x2_lb[i]) || x2 >= $(x2_ub[i]),
        )
    end
end

# Add rectangular obstacles to the model
add_rectangular_obstacles(model, model.variables["x1"], model.variables["x2"], x1_lb, x1_ub, x2_lb, x2_ub, N)

# Objective: Minimize time to reach the target position
CoreControls.@minimize(model, :(sum(1 for t in 0:$(N-1))))

# Define the algorithm (e.g., uniform grid)
algorithm = CoreControls.@uniformgrid(model)

# Solve the problem
solution = CoreControls.solve(model, algorithm; horizon=10, Δt=0.1)

# Print the solution (Optional)
##print(solution)

# Visualize the path (Optional)
#TODO: Add visualization
##Visualize(solution, plot="path", labels=["x1", "x2"])
