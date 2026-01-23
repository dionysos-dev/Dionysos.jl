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

include("../problems/simple_pendulum.jl");

concrete_problem = SimplePendulum.safety_problem(; objective = "safety_up") # safety_up, safety_down

hx = SVector(3*(pi/180.0), 0.05)

u0 = SVector(0.0)
hu = SVector(0.3)

periodic_dims = SVector(1)
periods = SVector(2*pi)
periodic_start = SVector(-pi)

tstep = 0.1

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
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("use_periodic_domain"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)

MOI.optimize!(optimizer);

# Get the results
abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"));
abstract_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem"));
abstract_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_controller"));
concrete_controller = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))
concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"));
abstraction_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")
total_time = MOI.get(optimizer, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

nstep = 100

x0 = SVector(UT.sample(concrete_problem.initial_set)...)
x_traj, u_traj = ST.get_closed_loop_trajectory(
    MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system")),
    concrete_controller,
    x0,
    nstep;
    wrap = ST.get_periodic_wrapper(periodic_dims, periods; start = periodic_start),
)

# Here we display the coordinate projection on the two first components of the state space along the trajectory.
fig = plot(; aspect_ratio = :equal, legend = false);
plot!(
    UT.set_in_period(concrete_problem.system.X, periodic_dims, periods, periodic_start);
    color = :grey,
    opacity = 1.0,
    label = "",
);
plot!(
    UT.set_in_period(concrete_problem.safe_set, periodic_dims, periods, periodic_start);
    color = :red,
    opacity = 0.4,
    label = "Safe set",
);
plot!(
    UT.set_in_period(concrete_problem.initial_set, periodic_dims, periods, periodic_start);
    color = :green,
    opacity = 0.8,
    label = "Initial set",
);
plot!(x_traj; ms = 2.0, arrows = false)
display(fig)

# ### For Visualization

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
