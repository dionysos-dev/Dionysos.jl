using StaticArrays, Plots, Printf
import HybridSystems as HS
using JuMP

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/FlowShopScheduling/flow_shop_scheduling_3D.jl")

# -----------------------------------------------------
# Problem
# -----------------------------------------------------

concrete_problem = FlowShopScheduling3D.problem()
concrete_system = concrete_problem.system

# -----------------------------------------------------
# Sub-subsolvers parameters
# -----------------------------------------------------

# Discretization parameters (dx, du, dt) for each task
discretization_parameters = [(0.5, 0.5, 0.2), (0.5, 0.5, 0.2), (0.5, 0.5, 0.2)]

# Create optimizer factories for each mode using UniformGridAbstraction
optimizer_list = [
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
]

# Create state and input grids for each mode
state_grid_1 = MP.GridFree(
    SVector(0.0, 0.0, 0.0),
    SVector(
        discretization_parameters[1][1],
        discretization_parameters[1][1],
        discretization_parameters[1][1],
    ),
)
input_grid_1 = MP.GridFree(
    SVector(0.0, 0.0, 0.0),
    SVector(
        discretization_parameters[1][2],
        discretization_parameters[1][2],
        discretization_parameters[1][2],
    ),
)
state_grid_2 = MP.GridFree(
    SVector(0.0, 0.0, 0.0),
    SVector(
        discretization_parameters[2][1],
        discretization_parameters[2][1],
        discretization_parameters[2][1],
    ),
)
input_grid_2 = MP.GridFree(
    SVector(0.0, 0.0, 0.0),
    SVector(
        discretization_parameters[2][2],
        discretization_parameters[2][2],
        discretization_parameters[2][2],
    ),
)
state_grid_3 = MP.GridFree(
    SVector(0.0, 0.0, 0.0),
    SVector(
        discretization_parameters[3][1],
        discretization_parameters[3][1],
        discretization_parameters[3][1],
    ),
)
input_grid_3 = MP.GridFree(
    SVector(0.0, 0.0, 0.0),
    SVector(
        discretization_parameters[3][2],
        discretization_parameters[3][2],
        discretization_parameters[3][2],
    ),
)

# Jacobian bounds for each mode (used in abstraction growth)
jacobian_bounds = FlowShopScheduling3D.jacobian_bounds()

print_level = 2

# Create optimizer parameters dictionary
optimizer_kwargs_dict = [
    Dict(
        "state_grid" => state_grid_1,
        "input_grid" => input_grid_1,
        "time_step" => discretization_parameters[1][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[1],
        "print_level" => print_level,
    ),
    Dict(
        "state_grid" => state_grid_2,
        "input_grid" => input_grid_2,
        "time_step" => discretization_parameters[2][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[2],
        "print_level" => print_level,
    ),
    Dict(
        "state_grid" => state_grid_3,
        "input_grid" => input_grid_3,
        "time_step" => discretization_parameters[3][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[3],
        "print_level" => print_level,
    ),
]

# -----------------------------------------------------
# Solver parameters
# -----------------------------------------------------

optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
    optimizer_kwargs_dict,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

# -----------------------------------------------------
# Solve
# -----------------------------------------------------

MOI.optimize!(optimizer)

# -----------------------------------------------------
# Get Results
# -----------------------------------------------------

abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

# -----------------------------------------------------
# Closed-loop Trajectory
# -----------------------------------------------------

reached(aug_x) = AB.HybridSystemAbstraction.reached(concrete_problem, aug_x)

initial_state = concrete_problem.initial_set
max_steps = 10000
tsteps = [0.2, 0.2, 0.2] # per mode
aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
    concrete_system,
    concrete_controller,
    tsteps,
    initial_state,
    max_steps;
    stopping = reached,
)

for (idx, (t, u)) in enumerate(zip(aug_x_traj, u_traj))
    println("[", idx, "] state: ", t, " - control applied: ", u)
end
println("Final state: ", aug_x_traj[end])
