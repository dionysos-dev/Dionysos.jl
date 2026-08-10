# ============================================================================
# Thermostat as a HYBRID system WITHOUT the time lift.
#
# The modes are plain physical systems (no time subsystem), so the symbolic
# abstraction carries no time axis: its augmented state is `(temperature, mode)`
# instead of `(temperature, time, mode)`. This yields a strictly smaller
# abstraction and a time-independent specification (reach a temperature band in
# either mode, with no deadline). Compare with `thermostat_hybrid_time_system.jl`,
# which keeps the time lift and a time-windowed target.
# ============================================================================

using StaticArrays
using Plots
import LazySets

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
        "thermostat_hybrid_system.jl",
    ),
)

# ------------------------------------------------------------
# 1) Problem (time-free: no clock, no time window)
# ------------------------------------------------------------

params = ThermostatHybridSystem.Params(; Ta = 18.0, alpha = 0.2, beta = 5.0)

concrete_problem = ThermostatHybridSystem.problem(;
    params = params,
    initial_temperature = 18.0,
    initial_mode = 1,                       # start OFF and cold; must switch ON to heat
    target = LazySets.Hyperrectangle(; low = SVector(21.0), high = SVector(23.0)),
)
concrete_system = concrete_problem.system

# ------------------------------------------------------------
# 2) Per-mode subsolver parameters (mode 1 = OFF, mode 2 = ON)
# ------------------------------------------------------------

Δt = 0.1        # integration step for the dynamics (NOT a clock — there is no time axis)
h = 0.05        # temperature grid step
η = 0.1         # input grid step

optimizer_list = [
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
    () -> AB.UniformGridAbstraction.Optimizer{Float64}(),
]

state_grid = MP.GridFree(SVector(0.0), SVector(h))
input_grid = MP.GridFree(SVector(0.0), SVector(η))
jacobian_bounds = ThermostatHybridSystem.jacobian_bounds(params)

make_kwargs(jac) = Dict(
    "state_grid" => state_grid,
    "input_grid" => input_grid,
    "time_step" => Δt,
    "approx_mode" => AB.UniformGridAbstraction.GROWTH,
    "jacobian_bound" => jac,
    "print_level" => 1,
)
optimizer_kwargs_dict = [make_kwargs(jacobian_bounds[1]), make_kwargs(jacobian_bounds[2])]

# ------------------------------------------------------------
# 3) Solve
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

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

# The abstraction is genuinely time-free: no mode is clock-lifted, states are (T, mode).
@assert all(m -> m isa DI.Symbolic.SymbolicModelList, abstract_system.mode_models)
println("Abstraction states (T, mode): ", DI.Symbolic.get_n_state(abstract_system))

# ------------------------------------------------------------
# 4) Closed-loop trajectory  (augmented state is `(T, mode)`)
# ------------------------------------------------------------

reached(aug_x) = AB.HybridSystemAbstraction.reached(concrete_problem, aug_x)

initial_state = concrete_problem.initial_set
tsteps = [Δt, Δt]
max_steps = 1000

aug_x_traj, u_traj = AB.HybridSystemAbstraction.get_closed_loop_trajectory(
    concrete_system,
    concrete_controller,
    tsteps,
    initial_state,
    max_steps;
    stopping = reached,
)

println("Reached target in $(length(aug_x_traj) - 1) steps; final state: ", aug_x_traj[end])

# ------------------------------------------------------------
# 5) Dashboard animation: thermometer (left) + mode / temperature / input over time
# ------------------------------------------------------------

# Switch inputs are "SWITCH ..." labels; map them to NaN so the numeric input panel skips them.
u_numeric = [u isa AbstractString ? SVector(NaN) : SVector{1}(Float64.(u)) for u in u_traj]

# Decompose the augmented `(T, mode)` states into a channelled trajectory (states = temperature,
# modes = OFF/ON); there is no clock, so the dashboard derives elapsed time from `Δt`.
trajectory = AB.HybridSystemAbstraction.channelled_trajectory(aug_x_traj, u_numeric)

system_plot! = ThermostatHybridSystem.system_plot!(; problem = concrete_problem)
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    trajectory;
    xdims = (1,),
    udims = (1,),
    Δt = Δt,
    fps = 10,
    ylabel_mode = "mode",
    title = "Time-free hybrid thermostat",
    ylabel_state = "T [°C]",
    ylabel_input = "heater u",
    # filename = "thermostat_hybrid_system.gif",
)
