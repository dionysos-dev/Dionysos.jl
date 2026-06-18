using StaticArrays
using MathematicalSystems
using Dionysos
using JuMP
import MathOptInterface as MOI
import LinearAlgebra as LA

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction
const OPDS = OP.DiscreteSystems

# ------------------------------------------------------------
# 1) Define a simple 2D continuous-time system: x' = u
# ------------------------------------------------------------

include("../problems/toy_problem.jl")

_X_ = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(4.0, 4.0))
_U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))

concrete_system = ToyProblem.system(; _X_ = _X_, _U_ = _U_)
jacobian_bound = ToyProblem.jacobian_bound()

# ------------------------------------------------------------
# 2) Build abstraction
# ------------------------------------------------------------

alternating_simulation_problem =
    DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

x0_grid = SVector(-2.0, -2.0)
hx = SVector(0.2, 0.2)
state_grid = MP.GridFree(x0_grid, hx)

u0_grid = SVector(-1.0, -1.0)
hu = SVector(0.5, 0.5)
input_grid = MP.GridFree(u0_grid, hu)

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    alternating_simulation_problem,
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), 0.3)

MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH, # CENTER_SIMULATION GROWTH
)

MOI.set(optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
MOI.set(optimizer, MOI.RawOptimizerAttribute("n_samples"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.Silent(), true)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

println("Abstraction built.")

# ------------------------------------------------------------
# 3) Define optimal control problem
# ------------------------------------------------------------

_I_ = UT.HyperRectangle(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

g11 = UT.HyperRectangle(SVector(-1.0, 3.0), SVector(-0.3, 3.7))
g12 = UT.HyperRectangle(SVector(1.0, 2.0), SVector(3.0, 3.7))
target_set = UT.LazySetUnion([g11, g12])

state_cost = nothing
transition_cost = nothing
time_horizon = 20

concrete_problem = DI.Problem.OptimalControlProblem(
    concrete_system,
    _I_,
    target_set,
    state_cost,
    transition_cost,
    time_horizon,
)

# ------------------------------------------------------------
# 4) Generate closed-loop trajectory
# ------------------------------------------------------------

trajectory_generator = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
    optimizer;
    initial_state = SVector(-1.65, -1.65),
    concrete = false,
    nstep = 50,
)

AB.set_problem!(trajectory_generator, concrete_problem)
AB.generate!(trajectory_generator)

closed_loop_traj = AB.get_trajectory(trajectory_generator)
success = AB.get_success(trajectory_generator)
solve_time = AB.get_solve_time(trajectory_generator)

println("Trajectory generated.")
println("Success: ", success)
println("Solve time: ", solve_time)

xtraj = closed_loop_traj.x
utraj = closed_loop_traj.u

# ------------------------------------------------------------
# 5) Improve/refine trajectory with MPPI
# ------------------------------------------------------------

Δt = 0.3
discrete_system = ST.discretize_continuous_system(
    concrete_problem.system, Δt
)
discrete_problem = DI.Problem.discretize_problem(
    concrete_problem,
    Δt
)

noise_sampler = function (rng, u, k)
    σ = 0.3
    return σ * randn(rng, length(u))
end

project_input = function (u)
    return SVector(
        clamp(u[1], -1.0, 1.0),
        clamp(u[2], -1.0, 1.0),
    )
end

trajectory_cost = function(problem, traj)
    xs = traj.x.seq
    us = traj.u.seq

    # Distance to the closest target component at each time.
    target_distances = [
        minimum(LA.norm(x - UT.get_center(g)) for g in problem.target_set.sets)
        for x in xs
    ]

    best_target_distance = minimum(target_distances)

    # Prefer reaching the target early.
    hit_idx = findfirst(x -> x ∈ problem.target_set, xs)
    target_bonus = hit_idx === nothing ? 0.0 : -1000.0 / hit_idx

    # Penalize leaving the domain.
    domain_violation_cost = sum(x -> x ∈ problem.system.X ? 0.0 : 1000.0, xs)

    # Small control regularization.
    control_cost = sum(LA.norm(u)^2 for u in us)

    return 100.0 * best_target_distance +
           0.01 * control_cost +
           domain_violation_cost +
           target_bonus
end

using Random
rng = Random.default_rng()

mppi_generator = AB.MPPITrajectoryGenerator.TrajectoryGenerator(
    ;
    rng = rng,
    seed_trajectory = closed_loop_traj,
    nstep = 50,
    nsamples = 1000,
    niter = 20,
    λ = 1.0,
    noise_sampler = noise_sampler,
    project_input = project_input,
    trajectory_cost = trajectory_cost,
    hard_constraint = false,
)

AB.set_problem!(mppi_generator, discrete_problem)
AB.generate!(mppi_generator)

mppi_traj = AB.get_trajectory(mppi_generator)
xtraj_mppi = mppi_traj.x

println("MPPI success: ", AB.get_success(mppi_generator))
println("MPPI solve time: ", AB.get_solve_time(mppi_generator))
println("MPPI diagnostics: ", AB.MPPITrajectoryGenerator.get_diagnostics(mppi_generator))

# ------------------------------------------------------------
# 6) Composite trajectory generator: use the closed-loop trajectory as a seed for MPPI
# ------------------------------------------------------------

combo_gen = AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
    trajectory_generator,
    mppi_generator;
    tstep = 0.3,
    num_substeps = 5,
)

AB.set_problem!(combo_gen, concrete_problem)
AB.generate!(combo_gen)

composite_traj = AB.get_trajectory(combo_gen)
composite_seed = AB.CompositeTrajectoryGenerator.get_seed(combo_gen)

# ------------------------------------------------------------
# 7) Plot
# ------------------------------------------------------------

using Plots

fig = plot(; aspect_ratio = :equal)

plot!(
    fig,
    concrete_problem;
    aspect_ratio = :equal,
)

plot!(
    fig,
    xtraj;
    color = :blue,
    dims = [1, 2],
    label = "Closed-loop trajectory",
)
plot!(
    fig,
    xtraj_mppi;
    color = :red,
    dims = [1, 2],
    label = "MPPI trajectory",
)

plot!(
    fig,
    composite_seed.x;
    color = :yellow,
    dims = [1, 2],
    label = "Composite seed trajectory",
)


plot!(
    fig,
    composite_traj.x;
    color = :green,
    dims = [1, 2],
    label = "Composite trajectory",
)

display(fig)