# Define state variables for the position and orientation
x1 = Variable(name="x1", domain=Reals())  # x position
x2 = Variable(name="x2", domain=Reals())  # y position
x3 = Variable(name="x3", domain=Reals())  # orientation

# Define control inputs with bounds representing U = [-1,1] × [-1,1]
u1 = Variable(name="u1", domain=Reals(), bounds=(-1, 1))  # rear wheel velocity
u2 = Variable(name="u2", domain=Reals(), bounds=(-1, 1))  # steering angle

# Define alpha as a function of u2
alpha_expr = :(atan(tan(u2) / 2))

# Dynamics of the system as a constraint
f_x1 = :(u1 * cos($alpha_expr + x3) * sec($alpha_expr))
f_x2 = :(u1 * sin($alpha_expr + x3) * sec($alpha_expr))
f_x3 = :(u1 * tan(u2))

# Add constraints for the system dynamics
Constraint(:($f_x1 == diff(x1)))  # Representing \( \dot{x_1} = f_x1 \)
Constraint(:($f_x2 == diff(x2)))  # Representing \( \dot{x_2} = f_x2 \)
Constraint(:($f_x3 == diff(x3)))  # Representing \( \dot{x_3} = f_x3 \)

# Define the objective function
# Assuming the objective is to minimize the time to reach the target
# The objective could be more complex depending on maze structure, but here we use a simplified time minimization example.
Objective(Minimize(), :(norm([x1 - target_x, x2 - target_y])))

# Capture solution and visualize state trajectory
solution = CaptureSolution(format="table", horizon=20)  # for a 20-step prediction or tracking horizon
Visualize(solution, plot="state_trajectory", labels=["x1", "x2", "x3"])
