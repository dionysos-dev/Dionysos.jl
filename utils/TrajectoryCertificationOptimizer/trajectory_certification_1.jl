using StaticArrays
using MathematicalSystems
using Dionysos
using JuMP
using Plots

import MathOptInterface as MOI
import LinearAlgebra as LA

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/toy_problem.jl")

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

_X_ = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(4.0, 4.0))
_U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))

concrete_system = ToyProblem.system(; _X_ = _X_, _U_ = _U_)
jacobian_bound = ToyProblem.jacobian_bound()

_I_ = UT.HyperRectangle(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

g11 = UT.HyperRectangle(SVector(-1.0, 3.0), SVector(-0.3, 3.7))
g12 = UT.HyperRectangle(SVector(1.0, 2.0), SVector(3.0, 3.7))
target_set = UT.set_union([g11, g12])

concrete_problem =
    DI.Problem.OptimalControlProblem(concrete_system, _I_, target_set, nothing, nothing, 20)

# ------------------------------------------------------------
# 2) Trajectory generator block 1: abstraction-based generator
# ------------------------------------------------------------

Δt = 0.3

alternating_simulation_problem =
    DI.Problem.AlternatingSimulationProblem(concrete_system, concrete_system.X)

state_grid = MP.GridFree(SVector(-2.0, -2.0), SVector(0.2, 0.2))
input_grid = MP.GridFree(SVector(-1.0, -1.0), SVector(0.5, 0.5))

abstraction_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    alternating_simulation_problem,
)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)
MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(abstraction_optimizer, MOI.Silent(), true)

trajectory_generator = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
    abstraction_optimizer;
    initial_state = SVector(-1.65, -1.65),
    concrete = false,
    nstep = 50,
)

# ------------------------------------------------------------
# 3) Trajectory generator block 2: MPPI generator
# ------------------------------------------------------------

noise_sampler = function (rng, u, k)
    σ = 0.3
    return σ * randn(rng, length(u))
end

project_input = function (u)
    return SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0))
end

trajectory_cost = function (problem, traj)
    xs = traj.x.seq
    us = traj.u.seq

    target_distances = [
        minimum(LA.norm(x - UT.get_center(g)) for g in problem.target_set.sets) for x in xs
    ]

    best_target_distance = minimum(target_distances)
    hit_idx = findfirst(x -> x ∈ problem.target_set, xs)
    target_bonus = hit_idx === nothing ? 0.0 : -1000.0 / hit_idx
    domain_violation_cost = sum(x -> x ∈ problem.system.X ? 0.0 : 1000.0, xs)
    control_cost = sum(LA.norm(u)^2 for u in us)

    return 100.0 * best_target_distance +
           0.01 * control_cost +
           domain_violation_cost +
           target_bonus
end

mppi_generator = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.default_rng(),
    seed_trajectory = nothing,
    nstep = 50,
    nsamples = 1000,
    niter = 20,
    λ = 1.0,
    noise_sampler = noise_sampler,
    project_input = project_input,
    trajectory_cost = trajectory_cost,
    hard_constraint = false,
)

# ------------------------------------------------------------
# 4) Composite trajectory generator
# ------------------------------------------------------------

combo_gen = AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
    trajectory_generator,
    mppi_generator;
    Δt = Δt,
    num_substeps = 5,
)

# ------------------------------------------------------------
# 5) Uniform-grid tube certifier
# ------------------------------------------------------------

const UGTC = AB.UniformGridTrajectoryCertifier

local_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(local_optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
MOI.set(local_optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(local_optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)

MOI.set(
    local_optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)

MOI.set(local_optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
MOI.set(local_optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(local_optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(local_optimizer, MOI.Silent(), true)

ug_certifier = UGTC.TrajectoryCertifier(;
    optimizer = local_optimizer,
    radius = 0.6,
    margin = 0.0,
    incl_mode = MP.INNER,
    n_between = 0,
    max_step = nothing,
    enforce_safe_max_step = true,
    handle_system_domain = true,
)

# ------------------------------------------------------------
# 6) Modular trajectory-generation + uniform-grid certification optimizer
# ------------------------------------------------------------

tc_optimizer = AB.TrajectoryCertificationOptimizer.Optimizer(combo_gen, ug_certifier)

MOI.set(tc_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

MOI.optimize!(tc_optimizer)

composite_traj = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("trajectory"))
controller = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("controller"))
success = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("success"))

println("Trajectory generation + uniform-grid certification success: ", success)
println(
    "Trajectory generation + uniform-grid certification time: ",
    MOI.get(tc_optimizer, MOI.SolveTimeSec()),
)

# For plotting
ug_certifier_optimizer = tc_optimizer.trajectory_certifier.optimizer
abstract_system =
    MOI.get(ug_certifier_optimizer, MOI.RawOptimizerAttribute("abstract_system"))
controllable_set =
    MOI.get(ug_certifier_optimizer, MOI.RawOptimizerAttribute("controllable_set"))
uncontrollable_set =
    MOI.get(ug_certifier_optimizer, MOI.RawOptimizerAttribute("uncontrollable_set"))
XMapping = SY.get_state_mapping(abstract_system)

# ------------------------------------------------------------
# 7) Plot
# ------------------------------------------------------------

fig = plot(; aspect_ratio = :equal)

plot!(fig, concrete_problem; aspect_ratio = :equal)

if composite_traj !== nothing
    X_tube = UGTC.build_tube(
        composite_traj.x,
        ug_certifier.radius;
        margin = ug_certifier.margin,
        n_between = ug_certifier.n_between,
        max_step = ug_certifier.max_step,
        enforce_safe_max_step = ug_certifier.enforce_safe_max_step,
        X_domain = ug_certifier.handle_system_domain ? concrete_problem.system.X : nothing,
    )

    plot!(fig, X_tube; alpha = 0.25, color = :blue, label = "Certified tube")

    plot!(
        (controllable_set, XMapping);
        color = :yellow,
        linecolor = :yellow,
        label = "Controllable set",
    )

    plot!(
        fig,
        composite_traj.x;
        color = :green,
        dims = [1, 2],
        linewidth = 2,
        label = "Generated trajectory",
    )
end

display(fig)
