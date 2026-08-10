using StaticArrays
using Dionysos
using Plots
using JuMP

import MathOptInterface as MOI
import LazySets

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
        "thermostat_system.jl",
    ),
)

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

params = ThermostatSystem.Params(;
    Ta = 20.0,
    alpha = 0.1, # disspativity coefficient: determine switching control frequency
    beta = 2.0,
)

concrete_problem = ThermostatSystem.reach_and_stay_problem(;
    params = params,
    _I_ = LazySets.Hyperrectangle(; low = SVector(18.0), high = SVector(18.5)),
    _T_ = LazySets.Hyperrectangle(; low = SVector(21.8), high = SVector(22.2)),
    _S_ = LazySets.Hyperrectangle(; low = SVector(17.0), high = SVector(25.0)),
)

concrete_system = concrete_problem.system

# ------------------------------------------------------------
# 2) Abstraction parameters
# ------------------------------------------------------------

Δt = 0.1

# Temperature grid size
h = SVector(0.05)

# Two modes:
#   u = 1 => heater ON
#   u = 2 => heater OFF
u0 = SVector(1)
hu = SVector(1)
input_grid = MP.GridFree(u0, hu)

jacobian_bound = ThermostatSystem.jacobian_bound(params)

# ------------------------------------------------------------
# 3) Build abstraction and synthesize controller
# ------------------------------------------------------------

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), h)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # GROWTH, CENTER_SIMULATION
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)

MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_update_interval"), 100)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.FastIndexedAutomatonList(n, m),
)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))

concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

println("success = ", success)

# ------------------------------------------------------------
# 4) Closed-loop simulation
# ------------------------------------------------------------

discrete_system = ST.discretize_continuous_system(concrete_system, Δt; num_substeps = 5)

x0 = SVector(18.2)
nstep = 50

traj = ST.get_closed_loop_trajectory(
    discrete_system,
    concrete_controller,
    x0,
    nstep;
    verbose = true,
)

# ------------------------------------------------------------
# 5) Dashboard animation
# ------------------------------------------------------------

system_plot! = ThermostatSystem.system_plot!(; problem = concrete_problem)

Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1,),
    udims = (1,),
    Δt = Δt,
    fps = 5,
    # filename = "thermostat_abstraction.mp4",
    title = "Thermostat abstraction controller",
    xlabel_state = "time [s]",
    ylabel_state = "temperature [°C]",
    xlabel_input = "time [s]",
    ylabel_input = "mode",
    ylims_state = (18.0, 24.0),
    ylims_input = (0.5, 2.5),
)
