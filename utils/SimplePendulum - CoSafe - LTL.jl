using StaticArrays, JuMP, Plots
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction

using Spot

# ------------------------------------------------------------
# 1) Define System
# ------------------------------------------------------------
include("../problems/simple_pendulum.jl");
_X_ = UT.HyperRectangle(SVector(-π, -5.5), SVector(π, 5.5))
_U_ = UT.LazySetMinus(
    UT.HyperRectangle(SVector(-4.5), SVector(4.5)),
    UT.HyperRectangle(SVector(-0.5), SVector(0.5)),
)
concrete_system = SimplePendulum.system(; l = 1.0, g = 9.81, _X_ = _X_, _U_ = _U_)

# ------------------------------------------------------------
# 2) Define co-safe LTL problem with LazySets-style labeling
# ------------------------------------------------------------

_I_ = UT.HyperRectangle(SVector(-5.0 * pi / 180.0, -0.2), SVector(5.0 * pi / 180.0, 0.2))

g1 = UT.HyperRectangle(
    SVector(pi - 10.0 * pi / 180.0, -1.0),
    SVector(pi + 15.0 * pi / 180.0, 1.0),
)

g2 = UT.HyperRectangle(
    SVector(pi/2.0-10.0 * pi / 180.0, -0.4),
    SVector(pi/2.0+10.0 * pi / 180.0, 0.4),
)

obs = UT.HyperRectangle(
    SVector(-pi + 16.0 * pi / 180.0, -5.5),
    SVector(-pi + 38.0 * pi / 180.0, 5.5),
)

φ = ltl"G(!obs) & F(g1 & F(g2))"

labeling = Dict{Symbol, Any}(:g1 => g1, :g2 => g2, :obs => obs)

ap_semantics = Dict{Symbol, Any}(:g1 => DO.INNER, :g2 => DO.INNER, :obs => DO.OUTER)

concrete_problem =
    PR.CoSafeLTLProblem(concrete_system, _I_, φ, labeling, ap_semantics, false)

# ------------------------------------------------------------
# 3) Define solver meta-parameters
# ------------------------------------------------------------

hx = SVector(3*(pi/180.0), 0.05)

u0 = SVector(0.0)
hu = SVector(0.3)

tstep = 0.1

periodic_dims = SVector(1)
periods = SVector(2*pi)
periodic_start = SVector(-pi)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), DO.GridFree(u0, hu))
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("jacobian_bound"),
    SimplePendulum.jacobian_bound(),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), tstep)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # GROWTH, CENTER_SIMULATION
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_domain"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("automaton_constructor"),
    (n, m) -> SY.NewIndexedAutomatonList(n, m),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
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
q0 = MOI.get(optimizer, MOI.RawOptimizerAttribute("qa0"))
nstep = 100

x_traj, u_traj, q_traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    q0,
    nstep;
    update_on_next = true,
    stopping = x -> false,
    wrap = ST.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# ------------------------------------------------------------
# 7) Plots
# ------------------------------------------------------------

using Plots
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
plot!(fig, x_traj; color = :blue, dims = [1, 2])
display(fig)

# ------------------------------------------------------------
# 8) Visualization
# ------------------------------------------------------------

using RigidBodyDynamics
using MeshCat, MeshCatMechanisms

# --- build mechanism/state from your URDF ---
urdf = joinpath(dirname(dirname(pathof(Dionysos))), "problems/pendulum/", "Pendulum.urdf")
mechanism = parse_urdf(urdf)
state = MechanismState(mechanism)
joint = first(joints(mechanism))

# --- build trajectory data ---
state_values = [x_traj.seq[i] for i in 1:ST.length(x_traj)]
ts = collect(0.0:tstep:((length(state_values) - 1) * tstep))

# --- visualizer ---
mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf))
vis = mvis.visualizer
open(vis)

# --- animation (no Interpolations) ---
fps = round(Int, 1 / tstep)
anim = MeshCat.Animation(vis; fps = fps)

for k in eachindex(ts)
    θ = state_values[k][1]
    set_configuration!(state, joint, θ)

    MeshCat.atframe(anim, k) do
        return MeshCatMechanisms.set_configuration!(mvis, configuration(state))
    end
end

MeshCat.setanimation!(vis, anim; play = true)
