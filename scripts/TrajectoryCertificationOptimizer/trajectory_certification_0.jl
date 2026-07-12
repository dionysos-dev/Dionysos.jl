import LazySets
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

include("../../problems/ToyProblem/toy_problem.jl")

_X_ = UT.box(SVector(-2.0, -2.0), SVector(4.0, 4.0))
_U_ = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))

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
MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(optimizer, MOI.Silent(), true)

MOI.optimize!(optimizer)

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"))
discrete_time_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("discrete_time_system"))

println("Abstraction built.")

# ------------------------------------------------------------
# 3) Define optimal control problem
# ------------------------------------------------------------

_I_ = UT.box(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

g11 = UT.box(SVector(-1.0, 3.0), SVector(-0.3, 3.7))
g12 = UT.box(SVector(1.0, 2.0), SVector(3.0, 3.7))
target_set = UT.set_union([g11, g12])

state_cost = nothing
trans_cost = nothing
time_horizon = 20

concrete_problem = DI.Problem.OptimalControlProblem(
    concrete_system,
    _I_,
    target_set,
    state_cost,
    trans_cost,
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
discrete_system = ST.discretize_continuous_system(concrete_problem.system, Δt)
discrete_problem = DI.Problem.discretize_problem(concrete_problem, Δt)

noise_sampler = function (rng, u, k)
    σ = 0.3
    return SVector(σ * randn(rng), σ * randn(rng))
end

project_input = function (u)
    return SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0))
end

trajectory_cost = function (problem, traj)
    xs = traj.x.seq
    us = traj.u.seq

    # Distance to the closest target component at each time.
    target_distances = [
        minimum(LA.norm(x - LazySets.center(g)) for g in problem.target_set.array)
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

mppi_generator = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
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
    Δt = 0.3,
    num_substeps = 5,
)

AB.set_problem!(combo_gen, concrete_problem)
AB.generate!(combo_gen)

composite_traj = AB.get_trajectory(combo_gen)
composite_seed = AB.CompositeTrajectoryGenerator.get_seed(combo_gen)

# ------------------------------------------------------------
# 7) Certify trajectory with ellipsoidal backward certifier
# ------------------------------------------------------------

function plot_ellipsoid_chain!(
    fig,
    cert_result;
    max_ellipsoids = 100,
    label_prefix = "Backward ellipsoid",
)
    cert_result === nothing && return fig
    cert_result.lmi_data === nothing && return fig
    !haskey(cert_result.lmi_data, :ellipsoids) && return fig

    ellipsoids = cert_result.lmi_data.ellipsoids
    isempty(ellipsoids) && return fig

    idxs = unique(
        round.(
            Int,
            range(1, length(ellipsoids); length = min(max_ellipsoids, length(ellipsoids))),
        ),
    )

    for (j, idx) in enumerate(idxs)
        plot!(
            fig,
            ellipsoids[idx];
            label = j == 1 ? label_prefix : "",
            linewidth = 1.5,
            linestyle = :dash,
            alpha = 0.7,
        )
    end

    return fig
end

using Symbolics

const EB = AB.EllipsoidalBackwardTrajectoryCertifier

Symbolics.@variables x[1:2]
Symbolics.@variables u[1:2]
Symbolics.@variables w[1:2]

Δt = 0.3

fsymbolic = [x[1] + Δt * (u[1] + w[1]), x[2] + Δt * (u[2] + w[2])]

Wformat = UT.box(SVector(0.0, 0.0), SVector(0.0, 0.0))

provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    collect(x),
    collect(u),
    collect(w),
    [0.0, 0.0],              # ΔW radius
    ST.format_input_set(_U_),
    ST.format_noise_set(Wformat),
)

adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    false,                  # enabled
    [0.2, 0.2],              # ΔX_initial
    [0.01, 0.01],            # ΔX_min
    [1.0, 1.0],              # ΔX_max
    [0.5, 0.5],              # ΔU_initial
    [0.01, 0.01],            # ΔU_min
    [1.0, 1.0],              # ΔU_max
    2.0,                    # growth
    1.05,                   # safety
    1,                      # max_iters
    1e-8,                   # atol
    false,                  # verbose
    [1.0],                  # search_scales
    :first_consistent,      # objective
    true,                   # keep_first_consistent
)

ellip_opts = EB.EllipsoidalBackwardOptions(;
    # These limit the size of the computed predecessor ellipsoid and controller deviation. Larger values make certification easier.
    maxδx = 30,
    maxδu = 1.0,
    λ = 0.05, # This weights the transition-cost objective versus ellipsoid volume. In your current formulation,
    terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.5^2,
    # These define the box used to compute the affine approximation and Lipschitz bound. Larger boxes are more conservative but safer; smaller boxes are less conservative but may not contain the ellipsoid/controller image.
    linearization_δx = [0.2, 0.2],
    linearization_δu = [0.5, 0.5],
    adaptive_boxes = adaptive_opts,
    use_log_det = false,
)

using Clarabel
backend = Clarabel.Optimizer

certifier = EB.TrajectoryCertifier(provider, backend, ellip_opts)

AB.set_problem!(certifier, concrete_problem)
AB.set_trajectory!(certifier, composite_traj)
AB.certify!(certifier)

println("Certification success: ", AB.get_success(certifier))
println("Certification time: ", AB.get_solve_time(certifier))

controller = AB.get_controller(certifier)
cert_result = EB.get_result(certifier)

# ------------------------------------------------------------
# 8) Trajectory-generation + certification optimizer
# ------------------------------------------------------------

certifier = EB.TrajectoryCertifier(provider, Clarabel.Optimizer, ellip_opts)

combo_gen = AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
    trajectory_generator,
    mppi_generator;
    Δt = 0.3,
    num_substeps = 5,
)

tc_optimizer = AB.TrajectoryCertificationOptimizer.Optimizer(combo_gen, certifier)

MOI.set(tc_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

MOI.optimize!(tc_optimizer)

composite_traj = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("trajectory"))
controller = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("controller"))

println(
    "Trajectory generation + certification success: ",
    MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("success")),
)

println(
    "Trajectory generation + certification time: ",
    MOI.get(tc_optimizer, MOI.SolveTimeSec()),
)

# ------------------------------------------------------------
# 8) Plot
# ------------------------------------------------------------

using Plots

fig = plot(; aspect_ratio = :equal)

plot!(fig, concrete_problem; aspect_ratio = :equal)

plot_ellipsoid_chain!(fig, cert_result; max_ellipsoids = 10)

plot!(fig, xtraj; color = :blue, dims = [1, 2], label = "Closed-loop trajectory")
plot!(fig, xtraj_mppi; color = :red, dims = [1, 2], label = "MPPI trajectory")

plot!(
    fig,
    composite_seed.x;
    color = :yellow,
    dims = [1, 2],
    label = "Composite seed trajectory",
)

plot!(fig, composite_traj.x; color = :green, dims = [1, 2], label = "Composite trajectory")

display(fig)
