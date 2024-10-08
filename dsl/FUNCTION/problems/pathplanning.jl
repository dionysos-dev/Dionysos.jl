# Define state variables for the positions (x, y) and orientation
X = Variables(names=["x1", "x2", "x3"], domain=Reals())

# Define control inputs with bounds representing U = [-1,1] × [-1,1] rear wheel velocity and steering angle
U = Variables(names=["u1", "u2"], domain=Reals(), bounds=(-1, 1))

# Dynamics of the system as a constraint
alpha_expr = :(atan(tan(U[2]) / 2))
f_x1 = :(U[1] * cos($alpha_expr + X[3]) * sec($alpha_expr))
f_x2 = :(U[1] * sin($alpha_expr + X[3]) * sec($alpha_expr))
f_x3 = :(U[1] * tan(U[2]))

#TODO: maybe define the a Dynamics struct to encapsulate the dynamics of the system directly instead of using the constraints
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
