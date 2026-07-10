using StaticArrays
using Plots
using Printf
import LazySets
using JuMP

import MathOptInterface as MOI
import HybridSystems as HS

using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/Thermostat/thermostat_hybrid_system.jl")

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

params = ThermostatHybridSystem.Params(; Ta = 18.0, alpha = 0.2, beta = 5.0)

concrete_problem = ThermostatHybridSystem.problem(;
    params = params,
    initial_temperature = 18.0,
    initial_mode = 1,
    target = UT.box(SVector(21.0), SVector(23.0)),
)

concrete_system = concrete_problem.system

# ------------------------------------------------------------
# 2) Subsolver parameters
# ------------------------------------------------------------

# Mode 1 = OFF
# Mode 2 = ON

Δt_off = 0.1
Δt_on = 0.1

h_off = 0.05
h_on = 0.05

η_off = 0.1
η_on = 0.1

optimizer_list = [
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
]

state_grid_off = MP.GridFree(SVector(0.0), SVector(h_off))
state_grid_on = MP.GridFree(SVector(0.0), SVector(h_on))

input_grid_off = MP.GridFree(SVector(0.0), SVector(η_off))
input_grid_on = MP.GridFree(SVector(0.0), SVector(η_on))

jacobian_bounds = ThermostatHybridSystem.jacobian_bounds(params)

print_level = 2

optimizer_kwargs_dict = [
    Dict(
        "state_grid" => state_grid_off,
        "input_grid" => input_grid_off,
        "time_step" => Δt_off,
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[1],
        "print_level" => print_level,
    ),
    Dict(
        "state_grid" => state_grid_on,
        "input_grid" => input_grid_on,
        "time_step" => Δt_on,
        "approx_mode" => AB.UniformGridAbstraction.GROWTH,
        "jacobian_bound" => jacobian_bounds[2],
        "print_level" => print_level,
    ),
]

# ------------------------------------------------------------
# 3) Hybrid abstraction optimizer
# ------------------------------------------------------------

optimizer = MOI.instantiate(AB.HybridSystemAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("optimizer_list"), optimizer_list)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("optimizer_kwargs_dict"),
    optimizer_kwargs_dict,
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 1)

# ------------------------------------------------------------
# 4) Solve
# ------------------------------------------------------------

MOI.optimize!(optimizer)

abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"))
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

# ------------------------------------------------------------
# 5) Closed-loop trajectory
# ------------------------------------------------------------

reached(aug_x) = AB.HybridSystemAbstraction.reached(concrete_problem, aug_x)

initial_state = concrete_problem.initial_set
max_steps = 1000

tsteps = [Δt_off, Δt_on]

aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
    concrete_system,
    concrete_controller,
    tsteps,
    initial_state,
    max_steps;
    stopping = reached,
)

println("Closed-loop trajectory:")
for (idx, (x, u)) in enumerate(zip(aug_x_traj, u_traj))
    println("[", idx, "] state: ", x, " - control applied: ", u)
end

println("Final state: ", aug_x_traj[end])

# ------------------------------------------------------------
# 6) Extract trajectory data
# ------------------------------------------------------------

temperatures = [aug_x[1][1] for aug_x in aug_x_traj]
times = [aug_x[2] for aug_x in aug_x_traj]
modes = [aug_x[3] for aug_x in aug_x_traj]

function control_value(u)
    if isa(u, String)
        return NaN
    elseif isa(u, Number)
        return Float64(u)
    elseif isa(u, AbstractVector) && length(u) >= 1
        return Float64(u[1])
    else
        return NaN
    end
end

controls = [control_value(u) for u in u_traj]
control_times = times[1:length(controls)]

target = concrete_problem.target_set[1][1]
target_lower = LazySets.low(target, 1)
target_upper = LazySets.high(target, 1)

switch_indices = [i for i in 2:length(modes) if modes[i] != modes[i - 1]]

# ------------------------------------------------------------
# 7) Plot temperature
# ------------------------------------------------------------

fig_temp = plot(
    times,
    temperatures;
    linewidth = 2,
    color = :black,
    label = "Temperature",
    title = "Hybrid thermostat: temperature",
    xlabel = "time [s]",
    ylabel = "T [°C]",
    legend = :outerright,
)

plot!(
    fig_temp,
    times,
    fill(target_lower, length(times));
    linestyle = :dash,
    linewidth = 2,
    label = "target lower",
)

plot!(
    fig_temp,
    times,
    fill(target_upper, length(times));
    linestyle = :dash,
    linewidth = 2,
    label = "target upper",
)

scatter!(
    fig_temp,
    [times[1]],
    [temperatures[1]];
    marker = :star5,
    markersize = 8,
    label = "initial",
)

if !isempty(switch_indices)
    scatter!(
        fig_temp,
        times[switch_indices],
        temperatures[switch_indices];
        marker = :diamond,
        markersize = 7,
        label = "switch",
    )
end

display(fig_temp)

# ------------------------------------------------------------
# 8) Plot hybrid mode
# ------------------------------------------------------------

fig_mode = plot(
    times,
    modes;
    seriestype = :steppost,
    linewidth = 2,
    label = "mode",
    title = "Hybrid thermostat: mode",
    xlabel = "time [s]",
    ylabel = "mode",
    yticks = ([1, 2], ["OFF", "ON"]),
    ylims = (0.5, 2.5),
    legend = :outerright,
)

display(fig_mode)

# ------------------------------------------------------------
# 9) Plot continuous input
# ------------------------------------------------------------

fig_input = plot(
    control_times,
    controls;
    seriestype = :steppost,
    linewidth = 2,
    label = "u",
    title = "Hybrid thermostat: continuous input",
    xlabel = "time [s]",
    ylabel = "heater power u",
    ylims = (-0.05, 1.05),
    legend = :outerright,
)

display(fig_input)

# ------------------------------------------------------------
# 10) Combined plot
# ------------------------------------------------------------

fig_all = plot(fig_temp, fig_mode, fig_input; layout = (3, 1), size = (1000, 900))

display(fig_all)

# savefig(fig_all, "thermostat_hybrid_system.png")

println("Saved figure to thermostat_hybrid_system.png")
