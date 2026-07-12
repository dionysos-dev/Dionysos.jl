using StaticArrays, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

using Spot

# ------------------------------------------------------------
# 1) Define System
# ------------------------------------------------------------
include("../../problems/Pendulum/simple_pendulum.jl");

_X_ = UT.box(SVector(-π, -5.5), SVector(π, 5.5))
_U_ = UT.set_minus(UT.box(SVector(-4.5), SVector(4.5)), UT.box(SVector(-0.5), SVector(0.5)))
concrete_system = SimplePendulum.system(; _X_ = _X_, _U_ = _U_)

# ------------------------------------------------------------
# 2) Define co-safe LTL problem with sets labeling
# ------------------------------------------------------------

_I_ = UT.box(SVector(-5.0 * pi / 180.0, -0.2), SVector(5.0 * pi / 180.0, 0.2))

g1 = UT.box(SVector(pi - 10.0 * pi / 180.0, -1.0), SVector(pi + 15.0 * pi / 180.0, 1.0))

g2 = UT.box(SVector(pi/2.0-10.0 * pi / 180.0, -0.4), SVector(pi/2.0+10.0 * pi / 180.0, 0.4))

obs = UT.box(SVector(-pi + 16.0 * pi / 180.0, -5.5), SVector(-pi + 38.0 * pi / 180.0, 5.5))

φ = ltl"G(!obs) & F(g1 & F(g2))"
spec = Dionysos.spot_stepper(φ)

labeling = Dict{Symbol, Any}(:g1 => g1, :g2 => g2, :obs => obs)

ap_semantics = Dict{Symbol, Any}(:g1 => MP.INNER, :g2 => MP.INNER, :obs => MP.OUTER)

concrete_problem = PR.CoSafeLTLProblem(concrete_system, _I_, spec, labeling, ap_semantics)

# ------------------------------------------------------------
# 3) Define solver meta-parameters
# ------------------------------------------------------------

hx = SVector(3*(pi/180.0), 0.05)

u0 = SVector(0.0)
hu = SVector(0.3)

Δt = 0.1

periodic_dims = SVector(1)
periods = SVector(2*pi)
periodic_start = SVector(-pi)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), MP.GridFree(u0, hu))
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    SimplePendulum.jacobian_bound(),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # GROWTH, CENTER_SIMULATION
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.FastIndexedAutomatonList(n, m),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

# ------------------------------------------------------------
# 4) Solve the problem
# ------------------------------------------------------------

MOI.optimize!(optimizer);

# ------------------------------------------------------------
# 5) Get the results
# ------------------------------------------------------------
success = MOI.get(optimizer, MOI.RawOptimizerAttribute("success"))
println("Co-safe LTL success: $success")

concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))

abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")
total_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

# ------------------------------------------------------------
# 6) Visualization
# ------------------------------------------------------------

x0 = SVector(UT.sample(concrete_problem.initial_set)...)
nstep = 100

traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    update_on_next = true,
    stopping = x -> false,
    wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# ------------------------------------------------------------
# 7) Plots
# ------------------------------------------------------------

φ_str = string(φ)
fig = plot(; aspect_ratio = :equal, title = "$φ_str")
concrete_system = concrete_problem.system
plot!(
    UT.set_in_period(concrete_system.X, periodic_dims, periods, periodic_start);
    color = :grey,
    hole_color = :black,
    opacity = 1.0,
    label = "",
);
plot!(
    concrete_problem;
    ap_colors = Dict(:g1 => :red, :g2 => :cyan, :obs => :black),
    aspect_ratio = :equal,
)
plot!(fig, traj; color = :blue, dims = [1, 2])
display(fig)

# ------------------------------------------------------------
# 8) Animation with dashboard
# ------------------------------------------------------------

system_plot! = SimplePendulum.system_plot!()
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    traj;
    xdims = (1, 2),      # phase plot θ vs ω
    udims = (1,),        # input over time
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
