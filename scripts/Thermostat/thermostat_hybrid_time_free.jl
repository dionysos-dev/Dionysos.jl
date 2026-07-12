# ============================================================================
# Thermostat as a HYBRID system WITHOUT the time lift.
#
# The modes are plain physical systems (no time subsystem), so the symbolic
# abstraction carries no time axis: its augmented state is `(temperature, mode)`
# instead of `(temperature, time, mode)`. This yields a strictly smaller
# abstraction and a time-independent specification (reach a temperature band in
# either mode, with no deadline). Compare with `thermostat_hybrid_system.jl`,
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

include("../../problems/Thermostat/thermostat_hybrid_system.jl")

# ------------------------------------------------------------
# 1) Problem (time-free: no clock, no time window)
# ------------------------------------------------------------

params = ThermostatHybridSystem.Params(; Ta = 18.0, alpha = 0.2, beta = 5.0)

concrete_problem = ThermostatHybridSystem.problem_time_free(;
    params = params,
    initial_temperature = 18.0,
    initial_mode = 1,                       # start OFF and cold; must switch ON to heat
    target = UT.box(SVector(21.0), SVector(23.0)),
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
# 5) Extract trajectory data  (no clock: reconstruct elapsed time from the steps)
# ------------------------------------------------------------

temperatures = [aug_x[1][1] for aug_x in aug_x_traj]
modes = [aug_x[2] for aug_x in aug_x_traj]     # mode is the 2nd (last) component, not the 3rd

elapsed = zeros(length(aug_x_traj))
for i in 2:length(aug_x_traj)
    elapsed[i] = elapsed[i - 1] + tsteps[modes[i - 1]]
end

control_value(u) =
    isa(u, String) ? NaN :
    isa(u, Number) ? Float64(u) : (isa(u, AbstractVector) ? Float64(u[1]) : NaN)
controls = [control_value(u) for u in u_traj]

target = concrete_problem.target_set.per_mode[1].set   # StateSpec over temperature
target_lower = LazySets.low(target, 1)
target_upper = LazySets.high(target, 1)
switch_indices = [i for i in 2:length(modes) if modes[i] != modes[i - 1]]

# ------------------------------------------------------------
# 6) Plots
# ------------------------------------------------------------

fig_temp = plot(
    elapsed,
    temperatures;
    linewidth = 2,
    color = :black,
    label = "temperature",
    title = "Time-free hybrid thermostat: temperature",
    xlabel = "elapsed time [s]",
    ylabel = "T [°C]",
    legend = :outerright,
)
plot!(
    fig_temp,
    elapsed,
    fill(target_lower, length(elapsed));
    linestyle = :dash,
    label = "target lower",
)
plot!(
    fig_temp,
    elapsed,
    fill(target_upper, length(elapsed));
    linestyle = :dash,
    label = "target upper",
)
scatter!(
    fig_temp,
    [elapsed[1]],
    [temperatures[1]];
    marker = :star5,
    markersize = 8,
    label = "initial",
)
if !isempty(switch_indices)
    scatter!(
        fig_temp,
        elapsed[switch_indices],
        temperatures[switch_indices];
        marker = :diamond,
        markersize = 7,
        label = "switch",
    )
end

fig_mode = plot(
    elapsed,
    modes;
    seriestype = :steppost,
    linewidth = 2,
    label = "mode",
    title = "Time-free hybrid thermostat: mode",
    xlabel = "elapsed time [s]",
    ylabel = "mode",
    yticks = ([1, 2], ["OFF", "ON"]),
    ylims = (0.5, 2.5),
    legend = :outerright,
)

fig_input = plot(
    elapsed[1:length(controls)],
    controls;
    seriestype = :steppost,
    linewidth = 2,
    label = "u",
    title = "Time-free hybrid thermostat: input",
    xlabel = "elapsed time [s]",
    ylabel = "heater power u",
    ylims = (-0.05, 1.05),
    legend = :outerright,
)

fig_all = plot(fig_temp, fig_mode, fig_input; layout = (3, 1), size = (1000, 900))
display(fig_all)
# savefig(fig_all, "thermostat_hybrid_time_free.png")
