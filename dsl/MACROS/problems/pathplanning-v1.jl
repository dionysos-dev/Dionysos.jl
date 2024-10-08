# Define Model with specified grid type
model = Control("MazeSolver")  # UniformGrid() by default

# Define state variables for the position and orientation
x1 = @variable(model, "x1", Reals())
x2 = @variable(model, "x2", Reals())
x3 = @variable(model, "x3", Reals())

# Define control inputs with bounds
u1 = @variable(model, "u1", Reals(), (-1, 1))
u2 = @variable(model, "u2", Reals(), (-1, 1))

# Define alpha as a function of u2
alpha_expr = :(atan(tan(u2) / 2))

# Dynamics of the system as constraints
f_x1 = :(u1 * cos($alpha_expr + x3) * sec($alpha_expr))
f_x2 = :(u1 * sin($alpha_expr + x3) * sec($alpha_expr))
f_x3 = :(u1 * tan(u2))


#TODO: Define the a Dynamics struct to encapsulate the dynamics of the system directly instead of using the constraints
# Add constraints for the system dynamics
@constraint(model, :($f_x1 == diff(x1)))
@constraint(model, :($f_x2 == diff(x2)))
@constraint(model, :($f_x3 == diff(x3)))

# Define the objective function
@objective(model, Minimize(), :(norm([x1 - target_x, x2 - target_y])))

# Solve the problem with the specified grid type
solution = model.solve(horizon=20)

# Visualize the state trajectory
Visualize(solution, plot="state_trajectory", labels=["x1", "x2", "x3"])
