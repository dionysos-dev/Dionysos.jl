# Define Model with specified grid type
model = Control("MazeSolver")  # UniformGrid() by default

# Define state variables for the position and orientation
X = @variable(model, ["x1", "x2", "x3"], Reals())

# Define control inputs with bounds
U = @variable(model, ["u1", "u2"], Reals(), (-1, 1))

# Define alpha as a function of u2
alpha_expr = :(atan(tan(u2) / 2))

# Dynamics of the system as constraints
f = [
    :(u1 * cos($alpha_expr + x3) * sec($alpha_expr)),
    :(u1 * sin($alpha_expr + x3) * sec($alpha_expr)),
    :(u1 * tan(u2))
]

#TODO: Define the a Dynamics struct to encapsulate the dynamics of the system directly instead of using the constraints
# Add constraints for the system dynamics
@constraint(model, f, diff) # Representing \( \dot{x} = f \) 

# Define the objective function
@objective(model, Minimize(), :(norm([X - Target])))

# Solve the problem with the specified grid type
solution = model.solve(horizon=20)

# Visualize the state trajectory
Visualize(solution, plot="state_trajectory", labels=["x1", "x2", "x3"])
