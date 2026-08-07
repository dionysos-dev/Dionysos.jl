using StaticArrays
using JuMP
using Plots

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const PR = DI.Problem

import MathOptInterface as MOI

include(
    joinpath(
        dirname(dirname(pathof(Dionysos))),
        "problems",
        "Pendulum",
        "simple_pendulum.jl",
    ),
)

# ------------------------------------------------------------
# Concrete reach-and-stay problem: swing up and stay near upright
# ------------------------------------------------------------

params = SimplePendulum.Params()

_X_ = UT.box(SVector(-π, -7.0), SVector(π, 7.0))
_U_ = UT.box(SVector(-3.5), SVector(3.5)) # SVector(-3.0), SVector(3.0) is with escaping
_I_ = UT.box(SVector(-5.0π / 180.0, -0.2), SVector(5.0π / 180.0, 0.2))
_T_ = UT.box(SVector(π - 15.0π / 180.0, -1.0), SVector(π + 15.0π / 180.0, 1.0))
_S_ = _X_

concrete_system = SimplePendulum.system(; params = params, _X_ = _X_, _U_ = _U_)

concrete_problem = PR.ReachAndStayProblem(concrete_system, _I_, _T_, _S_)

# ------------------------------------------------------------
# Abstraction parameters
# ------------------------------------------------------------

hx = SVector(3.0 * π / 180.0, 0.05)

u0 = SVector(0.0)
hu = SVector(0.3)

periodic_dims = SVector(1)
periods = SVector(2π)
periodic_start = SVector(-π)

Δt = 0.1

# ------------------------------------------------------------
# Uniform-grid abstraction + reach-and-stay synthesis
# ------------------------------------------------------------

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), MP.GridFree(u0, hu))

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    SimplePendulum.jacobian_bound(params),
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)

MOI.set(optimizer, MOI.RawOptimizerAttribute("early_stop"), true)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.FastIndexedAutomatonList(n, m),
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer)

# ------------------------------------------------------------
# Results
# ------------------------------------------------------------

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
winning_set = MOI.get(optimizer, MOI.RawOptimizerAttribute("winning_set"))

abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")

abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")

total_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

# ------------------------------------------------------------
# Closed-loop simulation
# ------------------------------------------------------------

target_set =
    UT.set_in_period(concrete_problem.target_set, periodic_dims, periods, periodic_start)

safe_set =
    UT.set_in_period(concrete_problem.safe_set, periodic_dims, periods, periodic_start)

nstep = 150
x0 = SVector(UT.sample(concrete_problem.initial_set)...)

traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# ------------------------------------------------------------
# Check finite simulated reach-and-stay behavior
# ------------------------------------------------------------

success_check = PR.trajectory_success(concrete_problem, traj)

println("Simulation eventually stays in target over sampled horizon: $(success_check)")

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
XMapping = SY.get_state_mapping(abstract_system)

fig = plot(; aspect_ratio = :equal)

plot!(
    UT.set_in_period(concrete_system.X, periodic_dims, periods, periodic_start);
    color = :grey,
    hole_color = :black,
    opacity = 1.0,
    label = "",
)

plot!((winning_set, XMapping); color = :blue, linecolor = :blue)

plot!(
    UT.set_in_period(concrete_problem.initial_set, periodic_dims, periods, periodic_start);
    color = :green,
    opacity = 0.2,
    label = "Initial set",
)

plot!(target_set; color = :red, opacity = 0.8, label = "Target set")

plot!(traj; ms = 2.0, label = "Closed-loop trajectory")

display(fig)

# ------------------------------------------------------------
# Animation with dashboard
# ------------------------------------------------------------

system_plot! = SimplePendulum.system_plot!(; params = params)

Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1, 2),
    udims = (1,),
    Δt = Δt,
    fps = 5,
    # filename = "simple_pendulum_dashboard.mp4",
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    xlabel_input = "time [s]",
    ylabel_input = "τ [Nm]",
    xlims_state = (-π, π),
    ylims_state = (-8, 8),
)
