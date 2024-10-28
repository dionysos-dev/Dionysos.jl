# traffic_light_control.jl

# Import necessary modules
include("../src/core/core.jl")
include("../src/core/solver.jl")

using .CoreControls
using .Solver
using JuMP
using Ipopt
using Plots

# Create a Control model
model = Control(
    name="TrafficLightControl",
    problem_type=Stabilization(),
    system_type=Auto()
)

# Parameters
@parameter(model, "t_green", 30)
@parameter(model, "t_yellow", 5)
@parameter(model, "t_red", 35)
@parameter(model, "a", 5)
@parameter(model, "d", 5)
@parameter(model, "Q_max", 100)
@parameter(model, "Q_desired", 50)
@parameter(model, "beta", 0.1)
@parameter(model, "N", 20)
@parameter(model, "Q_initial", 70)

# State Variable
@variable(model, "Q", StateVar(), Reals(), bounds=(0, Q_max))

# Mode Variable
@variable(model, "mode", ModeVar(), Integers(), bounds=(1, 3))

# Mode Constants
const MODE_GREEN = 1
const MODE_YELLOW = 2
const MODE_RED = 3

# Dynamics and Mode Transitions
for t in 0:N-1
    @constraint(model, :(mode_$t == MODE_GREEN ==> Q_$(t+1) == Q_$t - d))
    @constraint(model, :(mode_$t == MODE_YELLOW ==> Q_$(t+1) == Q_$t))
    @constraint(model, :(mode_$t == MODE_RED ==> Q_$(t+1) == Q_$t + a))
    # Mode transitions can be further defined based on timing
end

# Initial Condition
@constraint(model, :(Q_0 == Q_initial))
@constraint(model, :(mode_0 == MODE_RED))

# Objective Function
@objective(model, Minimize(), :(
    sum((Q_$t - Q_desired)^2 for t in 0:N) + beta * sum(mode_$t != mode_$(t-1) for t in 1:N)
))

# Print the model
print(model)

# Define algorithm (since hybrid, may need specialized solver)
algorithm = SampleBasedAlgorithm(num_samples=100)

# Solve
solution = model.solve(algorithm; horizon=N)

# Visualize the queue length and mode
time = solution.solution_data["time"]
Q_values = solution.solution_data["Q"]
mode_values = solution.solution_data["mode"]
plot(time, Q_values, label="Queue Length", xlabel="Time", ylabel="Queue", title="Queue Length Over Time")
gui()

plot(time, mode_values, label="Mode", xlabel="Time", ylabel="Mode", title="Traffic Light Mode Over Time")
gui()
