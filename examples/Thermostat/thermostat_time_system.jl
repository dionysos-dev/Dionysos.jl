using StaticArrays
using Dionysos
using Plots
using JuMP

import MathOptInterface as MOI

const DI = Dionysos
const ST = DI.System
const UT = DI.Utils
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "Thermostat",
        "thermostat_time_system.jl",
    ),
)

# ------------------------------------------------------------
# 1) Problem: reach the temperature band within a clock window [min_time, max_time],
#    over the time-lifted (temperature, time) continuous system. The controller heats up
#    and holds T in the band until the window opens at min_time.
#
#    A positive min_time needs a time grid finer than Δt (see the state_grid comment
#    below) — otherwise the abstraction's clock can spuriously "stall" and the window
#    becomes unreachable.
# ------------------------------------------------------------

params = ThermostatTimeSystem.Params(; Ta = 20.0, alpha = 0.1, beta = 2.0)

concrete_problem = ThermostatTimeSystem.problem(;
    params = params,
    _I_ = UT.box(SVector(18.0, 0.0), SVector(18.5, 0.05)),
    _T_target = UT.box(SVector(21.0), SVector(23.0)),
    min_time = 3.0,
    max_time = 4.0,
)
concrete_system = concrete_problem.system

# ------------------------------------------------------------
# 2) Abstraction parameters ((T, t) grid; the time step matches the clock rate ṫ = 1)
# ------------------------------------------------------------

Δt = 0.1

# Grid step in (T, t). The clock must advance by MORE than one time cell per step:
# the GROWTH over-approximation adds ~one cell of margin around the nominal `+Δt` shift,
# so if the t-cell width equalled Δt that margin would reach back into the source cell,
# creating a spurious "clock may not advance" self-loop that makes any timed (min_time > 0)
# target unreachable. With h_t = Δt/2 the clock jumps ~2 cells while the margin is ~1, so
# it always advances. The origin is offset by half a cell so the first time cell sits fully
# inside the domain rather than straddling the t = 0 edge.
h_t = Δt / 2
state_grid = MP.GridFree(SVector(0.0, h_t / 2), SVector(0.05, h_t))

u0 = SVector(1)                 # u = 1 => heater ON, u = 2 => heater OFF
hu = SVector(1)
input_grid = MP.GridFree(u0, hu)

jacobian_bound = ThermostatTimeSystem.jacobian_bound(params)

# ------------------------------------------------------------
# 3) Build abstraction and synthesize controller
# ------------------------------------------------------------

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_update_interval"), 1000)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.FastIndexedAutomatonList(n, m),
)

MOI.optimize!(optimizer)

success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
println("success = ", success)

# ------------------------------------------------------------
# 4) Closed-loop simulation  (state is (T, t); the clock advances by Δt each step)
# ------------------------------------------------------------

discrete_system = ST.discretize_continuous_system(concrete_system, Δt; num_substeps = 5)

x0 = SVector(18.2, 0.0)
nstep = 120

traj = ST.get_closed_loop_trajectory(
    discrete_system,
    concrete_controller,
    x0,
    nstep;
    stopping = x -> x ∈ concrete_problem.target_set,
    verbose = true,
)

println(
    "Reached target in $(length(ST.states(traj)) - 1) steps; final state: ",
    ST.states(traj)[end],
)

# ------------------------------------------------------------
# 5) Dashboard animation
# ------------------------------------------------------------

system_plot! = ThermostatTimeSystem.system_plot!(; problem = concrete_problem)
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1,),
    udims = (1,),
    Δt = Δt,
    fps = 5,
    title = "Time-lifted thermostat",
    ylabel_state = "T [°C]",
    ylabel_input = "mode",
    ylims_state = (17.0, 25.0),
    ylims_input = (0.5, 2.5),
    # filename = "thermostat_time_system.gif",
)
