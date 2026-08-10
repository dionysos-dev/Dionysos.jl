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

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "Thermostat",
        "thermostat_hybrid_time_system.jl",
    ),
)

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

params = ThermostatHybridTimeSystem.Params(; Ta = 18.0, alpha = 0.2, beta = 5.0)

concrete_problem = ThermostatHybridTimeSystem.problem(;
    params = params,
    initial_temperature = 18.0,
    initial_mode = 1,
    target = LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0)),
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

jacobian_bounds = ThermostatHybridTimeSystem.jacobian_bounds(params)

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

# ------------------------------------------------------------
# 6) Dashboard animation: thermometer (left) + mode / temperature / input over time
# ------------------------------------------------------------

# Switch inputs are "SWITCH ..." labels; map them to NaN so the numeric input panel skips them.
u_numeric = [u isa AbstractString ? SVector(NaN) : SVector{1}(Float64.(u)) for u in u_traj]

# Decompose the augmented `(T, time, mode)` states into a channelled trajectory (states =
# temperature, times = clock, modes = OFF/ON); the dashboard reads the channels directly.
trajectory = AB.HybridSystemAbstraction.channelled_trajectory(aug_x_traj, u_numeric)

system_plot! = ThermostatHybridTimeSystem.system_plot!(; problem = concrete_problem)
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    trajectory;
    xdims = (1,),
    udims = (1,),
    Δt = Δt_off,
    fps = 30,
    frame_step = 3, # render every 3rd step (~166 steps -> ~56 frames) to animate faster
    ylabel_mode = "mode",
    title = "Time-lifted hybrid thermostat",
    ylabel_state = "T [°C]",
    ylabel_input = "heater u",
    # filename = "thermostat_hybrid_time_system.gif",
)
