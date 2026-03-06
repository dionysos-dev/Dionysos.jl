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

include("../../problems/FlowShopScheduling/flow_shop_scheduling_1D.jl")

# -----------------------------------------------------
# Problem
# -----------------------------------------------------

concrete_problem = FlowShopScheduling1D.problem()
concrete_system = concrete_problem.system

# -----------------------------------------------------
# Sub-subsolvers parameters
# -----------------------------------------------------

# Discretization parameters (dx, du, dt) for each task
discretization_parameters = [
    (0.25, 0.25, 0.2),
    (0.25, 0.5, 0.1),
    (0.1, 0.1, 0.2),
    (0.1, 0.1, 0.1),
    (0.2, 0.25, 0.1),
]

# Create optimizer factories for each mode using UniformGridAbstraction
optimizer_list = [
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
]

# Create state and input grids for each mode
state_grid_1 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[1][1]))
input_grid_1 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[1][2]))
state_grid_2 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[2][1]))
input_grid_2 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[2][2]))
state_grid_3 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[3][1]))
input_grid_3 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[3][2]))
state_grid_4 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[4][1]))
input_grid_4 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[4][2]))
state_grid_5 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[5][1]))
input_grid_5 = MP.GridFree(SVector(0.0), SVector(discretization_parameters[5][2]))

# Jacobian bounds for each mode (used in abstraction growth)
jacobian_bounds = FlowShopScheduling1D.jacobian_bounds()

# Create optimizer parameters dictionary
optimizer_kwargs_dict = [
    Dict(
        "state_grid" => state_grid_1,
        "input_grid" => input_grid_1,
        "time_step" => discretization_parameters[1][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[1],
    ),
    Dict(
        "state_grid" => state_grid_2,
        "input_grid" => input_grid_2,
        "time_step" => discretization_parameters[2][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[2],
    ),
    Dict(
        "state_grid" => state_grid_3,
        "input_grid" => input_grid_3,
        "time_step" => discretization_parameters[3][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[3],
    ),
    Dict(
        "state_grid" => state_grid_4,
        "input_grid" => input_grid_4,
        "time_step" => discretization_parameters[4][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[4],
    ),
    Dict(
        "state_grid" => state_grid_5,
        "input_grid" => input_grid_5,
        "time_step" => discretization_parameters[5][3],
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[5],
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
tsteps = [0.2, 0.1, 0.2, 0.1, 0.1] # per mode
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

# -----------------------------------------------------
# Plot results
# -----------------------------------------------------

# Extract trajectory data
t_traj = [s[2] for s in aug_x_traj]
x_traj = [s[1][1] for s in aug_x_traj]
k_traj = [s[3] for s in aug_x_traj]

# Plot limits
t_min, t_max = minimum(t_traj), maximum(t_traj)
x_min, x_max = minimum(x_traj), maximum(x_traj)
t_margin = (t_max - t_min) * 0.05
x_margin = (x_max - x_min) * 0.1

fig = plot(;
    title = "System state over time (1D flowshop scheduling)",
    xlabel = "Time (s)",
    ylabel = "System state (x)",
    xlims = (t_min - t_margin, t_max + t_margin),
    ylims = (x_min - x_margin, x_max + x_margin),
    legend = :outerright,
    legendtitle = "Legend",
    size = (1400, 900),
    dpi = 150,
)

# Draw acceptance regions (guards)
region_colors = palette(:tab20)
auto = concrete_system.automaton
for (i, t) in enumerate(HS.transitions(auto))
    guard = HS.guard(concrete_system, t)
    if guard !== nothing
        x_accepted = [guard.lb[1], guard.ub[1]]
        t_accepted = [guard.lb[2], guard.ub[2]]
        color_task = region_colors[(i - 1) % length(region_colors) + 1]
        plot!(
            [t_accepted[1], t_accepted[2], t_accepted[2], t_accepted[1], t_accepted[1]],
            [x_accepted[1], x_accepted[1], x_accepted[2], x_accepted[2], x_accepted[1]];
            fill = (0, 0.15),
            fillcolor = color_task,
            alpha = 0.3,
            linecolor = color_task,
            linewidth = 2,
            linestyle = :dash,
            label = "Guard $(i): x∈[$(x_accepted[1]),$(x_accepted[2])], t∈[$(t_accepted[1]),$(t_accepted[2])]",
        )
        # Add grid for each task
        dx, du, dt = discretization_parameters[i]
        # Vertical grid lines
        t_grid = t_accepted[1]:dt:t_accepted[2]
        for t in t_grid
            if t >= t_min && t <= t_max
                vline!(
                    [t];
                    color = color_task,
                    alpha = 0.25,
                    linewidth = 1,
                    linestyle = :dot,
                    label = false,
                )
            end
        end
        # Horizontal grid lines
        x_grid = (floor(x_accepted[1] / dx) * dx):dx:(ceil(x_accepted[2] / dx) * dx)
        for x in x_grid
            if x >= x_min && x <= x_max
                hline!(
                    [x];
                    color = color_task,
                    alpha = 0.25,
                    linewidth = 1,
                    linestyle = :dot,
                    label = false,
                )
            end
        end
    end
end

# Plot the final target region
target_set = concrete_problem.target_set
Xs_target = target_set[1]
Ts_target = target_set[2]
final_x_target = Xs_target[1]
final_t_target = Ts_target[1]
color_final = :magenta
if !isnothing(final_x_target) && !isnothing(final_t_target)
    x_accepted = [final_x_target.lb[1], final_x_target.ub[1]]
    t_accepted = [final_t_target.lb[1], final_t_target.ub[1]]
    plot!(
        [t_accepted[1], t_accepted[2], t_accepted[2], t_accepted[1], t_accepted[1]],
        [x_accepted[1], x_accepted[1], x_accepted[2], x_accepted[2], x_accepted[1]];
        fill = (0, 0.25),
        fillcolor = color_final,
        alpha = 0.4,
        linecolor = color_final,
        linewidth = 3,
        linestyle = :solid,
        label = "Final target (problem target)",
    )
end

# Plot main trajectory
plot!(
    t_traj,
    x_traj;
    color = :black,
    linewidth = 2,
    label = "Trajectory",
    linestyle = :solid,
    alpha = 0.8,
)

## Initial point and switches
scatter!(
    [t_traj[1]],
    [x_traj[1]];
    color = :green,
    marker = :star5,
    markersize = 8,
    label = "Initial state",
)
switch_indices = [i for i in 2:length(k_traj) if k_traj[i] != k_traj[i - 1]]
if !isempty(switch_indices)
    t_switches = [t_traj[i] for i in switch_indices]
    x_switches = [x_traj[i] for i in switch_indices]
    scatter!(
        t_switches,
        x_switches;
        color = :red,
        marker = :diamond,
        markersize = 8,
        label = "Switch",
    )
    for i in switch_indices
        plot!(
            [t_traj[i - 1], t_traj[i]],
            [x_traj[i - 1], x_traj[i]];
            color = :red,
            linewidth = 3,
            linestyle = :dash,
            alpha = 0.9,
            label = false,
        )
    end
end

## Task starts
if !isempty(switch_indices)
    t_starts = [t_traj[i] for i in switch_indices]
    x_starts = [x_traj[i] for i in switch_indices]
    scatter!(
        t_starts,
        x_starts;
        color = :blue,
        marker = :circle,
        markersize = 5,
        label = "Task start",
    )
end

## Control annotations (sampled)
n_annotations = min(15, length(u_traj))
step = max(1, length(u_traj) ÷ n_annotations)
for i in 1:step:length(u_traj)
    if i <= length(t_traj) && u_traj[i] != "switch" && u_traj[i] != "SWITCH"
        if isa(u_traj[i], Number)
            control_str = Printf.@sprintf("u = %.2f", Float64(u_traj[i]))
        elseif isa(u_traj[i], Vector) && length(u_traj[i]) > 0
            control_values =
                join([Printf.@sprintf("%.2f", Float64(u)) for u in u_traj[i]], ", ")
            control_str = "u = $control_values"
        elseif isa(u_traj[i], AbstractArray) && length(u_traj[i]) == 1
            control_str = Printf.@sprintf("u = %.2f", Float64(u_traj[i][1]))
        else
            control_str = "u = $(u_traj[i])"
        end
        t_pos = t_traj[i]
        x_pos = x_traj[i] + x_margin * 0.3
        annotate!(t_pos, x_pos, text(control_str, :black, 8, :left))
    end
end

display(fig)
