# inventory_management.jl

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
    name="InventoryManagement",
    problem_type=Safety(),
    system_type=Discrete()
)

# Parameters
@parameter(model, "I_min", 50.0)
@parameter(model, "I_max", 500.0)
@parameter(model, "S_max", 100.0)
@parameter(model, "D_min", 20.0)
@parameter(model, "D_max", 80.0)
@parameter(model, "N", 10)
@parameter(model, "I_initial", 200.0)
@parameter(model, "holding_cost", 1.0)
@parameter(model, "ordering_cost", 2.0)

# State Variable
@variable(model, "I", StateVar(), Reals(), bounds=(I_min, I_max))

# Control Input
@variable(model, "S", InputVar(), Reals(), bounds=(0, S_max))

# Demand as a parameter (can be uncertain within bounds)
@parameter(model, "D", 50.0)  # For simplicity, using a fixed demand

# Dynamics
for k in 0:N-1
    @constraint(model, :(I_$(k+1) == I_$k + S_$k - D))
end

# Safety Stock Constraint
for k in 0:N
    @constraint(model, :(I_$k >= I_min))
end

# Initial Condition
@constraint(model, :(I_0 == I_initial))

# Objective Function
@objective(model, Minimize(), :(sum(holding_cost * I_$k + ordering_cost * S_$k for k in 0:$(N-1))))

# Print the model
print(model)

# Define algorithm (not specified, so we use default)
algorithm = UniformGridAlgorithm(grid_size=10)

# Solve
solution = model.solve(algorithm; horizon=N)

# Visualize the inventory levels
time = solution.solution_data["time"]
I_values = solution.solution_data["I"]
plot(time, I_values, label="Inventory Level", xlabel="Time", ylabel="Inventory", title="Inventory Levels Over Time")
#gui()
