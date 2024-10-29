# Define the thermostat control model with an added control input
model = @control(
    name="ThermostatControlSystemWithControl",
    problem_type=Safety(),
    system_type=Auto()
)

# Parameters
@parameter(model, "T_set", 22.0)          # Desired temperature in degrees Celsius
@parameter(model, "delta", 1.0)           # Temperature tolerance
@parameter(model, "k", 0.5)               # Maximum heating rate
@parameter(model, "alpha", 0.1)           # Cooling rate when in idle mode

# State Variable
@statevar(model, "T", Reals())            # Room temperature in real domain

# Control Variable (Heating Intensity)
@inputvar(model, "u", Reals(), bounds=(0, 1))  # Control input with bounds [0, 1]

# Mode Variable (Heating or Idle)
@modevar(model, "mode", 2)                # Two modes: Heating (1), Idle (0)

# Dynamics
@constraint(model, :(mode == 1 => { dot(T) == k * u } ))      # Heating mode dynamics with control input
@constraint(model, :(mode == 0 => { dot(T) == -alpha } ))     # Idle mode dynamics (ambient cooling)

# Safety Constraints
@constraint(model, T >= T_set - delta)            # Lower bound for temperature
@constraint(model, T <= T_set + delta)            # Upper bound for temperature

# define algorithm
algo = @uniformgrid(model)

# Solve
solution = solve(model, algo)
Visualize(solution)                               # Visualize solution for control switching and temperature changes
