using .CoreControls

# Example Usage
# Declare variables
x = Variable("x", Reals())
u = Variable("u", PositiveReals(), bounds=(0, 10))

# Define constraints
c1 = Constraint(:(x + u == 10))
c2 = Constraint(:(x - 2*u >= 0))

# Define an objective
objective = Objective(Minimize(), :(x^2 + u^2))

# Capture and visualize solution
solution = CaptureSolution(format="table", horizon=10)
Visualize(solution, plot="state_trajectory")