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

include("../../problems/CartPendulum/cart_pendulum.jl")

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

params =
    CartPendulum.Params(; M = 1.0, m = 0.2, l = 1.0, J = 0.0, c = 0.1, γ = 0.05, g = 9.81)

concrete_problem =
    CartPendulum.optimal_control_problem(; params = params, objective = "swing_up")

concrete_system = concrete_problem.system

# ------------------------------------------------------------
# 2) Abstraction parameters
# ------------------------------------------------------------

Δt = 0.1

periodic_dims = SVector(2)
periodic_periods = SVector(2π)
periodic_start = SVector(-π)

h = SVector(0.25, 10.0π / 180.0, 0.5, 0.5)

input_grid = MP.GridFree(SVector(0.0), SVector(1.0))

jacobian_bound = CartPendulum.jacobian_bound(params)

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
    AB.UniformGridAbstraction.CENTER_SIMULATION, # GROWTH, CENTER_SIMULATION
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)

MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periodic_periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)

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

# ------------------------------------------------------------
# 5) Closed-loop simulation
# ------------------------------------------------------------

discrete_system = ST.discretize_continuous_system(concrete_system, Δt; num_substeps = 5)

wrap = UT.get_periodic_wrapper(periodic_dims, periodic_periods; start = periodic_start)

reached(x) = wrap(x) ∈ concrete_problem.target_set

x0 = SVector(0.0, 0.0, 0.0, 0.0)
nstep = 500

traj = ST.get_closed_loop_trajectory(
    discrete_system,
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
    wrap = wrap,
)

# ------------------------------------------------------------
# 6) Plot θ-ω trajectory
# ------------------------------------------------------------

fig = plot(; aspect_ratio = :equal)

Xplot = UT.set_in_period(
    concrete_problem.system.X,
    periodic_dims,
    periodic_periods,
    periodic_start,
)

Iplot = UT.set_in_period(
    concrete_problem.initial_set,
    periodic_dims,
    periodic_periods,
    periodic_start,
)

Tplot = UT.set_in_period(
    concrete_problem.target_set,
    periodic_dims,
    periodic_periods,
    periodic_start,
)

plot!(fig, Xplot; dims = [2, 4], color = :grey, opacity = 0.15, label = "X")

plot!(fig, Iplot; dims = [2, 4], color = :green, opacity = 0.3, label = "Initial")

plot!(fig, Tplot; dims = [2, 4], color = :red, opacity = 0.3, label = "Target")

plot!(fig, traj; dims = [2, 4], linewidth = 2, label = "closed-loop θ-ω")

xlabel!(fig, "θ")
ylabel!(fig, "ω")

display(fig)

# ------------------------------------------------------------
# 7) Dashboard animation
# ------------------------------------------------------------

system_plot! =
    CartPendulum.system_plot!(; params = params, xlims = (-5.5, 5.5), ylims = (-1.4, 1.4))

Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (2, 4),
    udims = (1,),
    Δt = Δt,
    fps = 5,
    # filename = "cart_pendulum_abstraction.mp4",
    title = "Cart-pendulum abstraction controller",
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    xlabel_input = "time [s]",
    ylabel_input = "F [N]",
    xlims_state = (-π, π),
    ylims_state = (-5.0, 5.0),
    ylims_input = (-10.5, 10.5),
)
