# Define the thermostat control model
# main.jl

# Import necessary modules
include("../src/core/core.jl")

using .CoreControls


model = @control(
    name="ThermostatControlSystem",
    problem_type=Safety(),
    system_type=Hybrid()
)

# Parameters
@parameter(model, "T_set", 22.0)          # Desired temperature in degrees Celsius
@parameter(model, "delta", 1.0)           # Temperature tolerance
@parameter(model, "heating_rate", 0.5)    # Rate of temperature increase per time unit

# State Variable
@statevar(model, "T", Reals())            # Room temperature in real domain

# Mode Variable (Heating or Idle)
@modevar(model, "mode", 2)                # Two modes: Heating (1), Idle (0)

# Dynamics
@constraint(model, :(mode == 1 ==> dot(T) == heating_rate))     # Heating mode dynamics
@constraint(model, :(mode == 0 ==> dot(T) == 0))                 # Idle mode dynamics (or gradual cooling if needed)

# Safety Constraints
@constraint(model, T >= T_set - delta)    # Lower bound for temperature
@constraint(model, T <= T_set + delta)    # Upper bound for temperature

# Initial Conditions
T_start = 20.0
@constraint(model, :(T_0 == T_start))

# Terminal Conditions
T_goal = 22.0
@constraint(model, :(T_1 == T_goal))

# define algorithm
algo = @uniformgrid(model)

# Solve
solution = solve(model, algo)

Visualize(solution)                       # Visualize solution for control switching and temperature changes
