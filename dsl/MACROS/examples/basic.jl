# Example Usage
model = Control("Basic Example") 

# Declare variables
x = @variable(model, "x", Reals())
u = @variable(model, "u", PositiveReals(), bounds=(0, 10))

# Define constraints
@constraint(model, :(x + u == 10))
@constraint(model, :(x - 2*u >= 0))

# Define an objective
@objective(model, Minimize(), :(x^2 + u^2))

# Capture and visualize solution
solution = model.solve()

Visualize(solution, plot="state_trajectory")