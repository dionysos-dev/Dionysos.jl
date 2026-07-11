using StaticArrays
using MathematicalSystems
using Dionysos
import LazySets
using JuMP
using Random
using Symbolics
using Clarabel
using Plots

import MathOptInterface as MOI
import LinearAlgebra as LA

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

include("../../problems/toy_problem.jl")

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

_X_ = UT.box(SVector(-3.0, -3.0), SVector(4.0, 4.0))
_U_ = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))

concrete_system = ToyProblem.system(; _X_ = _X_, _U_ = _U_)
jacobian_bound = ToyProblem.jacobian_bound()

_I_ = UT.box(SVector(-1.7, -1.7), SVector(-1.6, -1.6))

g11 = UT.box(SVector(-1.0, 3.0), SVector(-0.3, 3.7))
g12 = UT.box(SVector(1.0, 2.0), SVector(3.0, 3.7))
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
    AB.UniformGridAbstraction.CENTER_SIMULATION, # GROWTH, CENTER_SIMULATION
)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)
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
    return SVector(σ * randn(rng), σ * randn(rng))
end

project_input = function (u)
    return SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0))
end

trajectory_cost = function (problem, traj)
    xs = traj.x.seq
    us = traj.u.seq

    target_distances = [
        minimum(LA.norm(x - LazySets.center(g)) for g in problem.target_set.array)
        for x in xs
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
# 5) Ellipsoidal backward certifier
# ------------------------------------------------------------

const EB = AB.EllipsoidalBackwardTrajectoryCertifier

Symbolics.@variables x[1:2]
Symbolics.@variables u[1:2]
Symbolics.@variables w[1:2]

# Defines symbolic state, input, and disturbance variables
fsymbolic = [x[1] + Δt * (u[1] + w[1]), x[2] + Δt * (u[2] + w[2])]

# Defines zero disturbance set
Wformat = UT.box(SVector(0.0, 0.0), SVector(0.0, 0.0))

# This object tells the certifier how to build local affine approximations along the trajectory.
provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,                   # symbolic dynamics
    collect(x),                  # symbolic state variables
    collect(u),                  # symbolic input variables
    collect(w),                  # symbolic disturbance variables
    [0.0, 0.0],                  # disturbance radius ΔW
    ST.format_input_set(_U_),    # input constraints for LMI
    ST.format_noise_set(Wformat), # disturbance vertices for LMI
)

# These is only used if the first argument is true
adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    false,        # enabled
    [0.2, 0.2],  # initial state linearization box radius
    [0.01, 0.01], # minimum state box radius
    [2.0, 2.0] * 10000,  # maximum state box radius
    [0.5, 0.5],  # initial input linearization box radius
    [0.01, 0.01], # minimum input box radius
    [2.0, 2.0],  # maximum input box radius
    2.0,         # growth factor if box is too small
    1.05,        # safety factor around required box size
    10,          # maximum adaptive iterations
    1e-8,        # numerical tolerance
    false,       # verbose
    [0.5, 0.75, 1.0, 1.25, 1.5], # candidate rescaling factors
    :first_consistent,
    true,        # accept first consistent box
)

trans_cost = UT.QuadraticStateControlFunction(
    zeros(2, 2),              # Q: state cost
    Matrix{Float64}(LA.I, 2, 2), # R: input cost
    zeros(2, 2),              # N
    zeros(2),                 # q
    zeros(2),                 # r
    0.0,                      # v
)

# Note: The toy system considered here is already affine:
# therefore the local affine approximation is exact. The Jacobians are constant, the Lipschitz remainder is identically zero, and the nonlinear error bounds vanish. Consequently, the parameters linearization_δx and linearization_δu do not affect the accuracy of the model approximation itself. Their only role is through the consistency checks used by the certification algorithm.
# In particular, when adaptive linearization boxes are enabled, the algorithm computes the state and input radii required to certify a transition. Since the dynamics are exactly linear, these radii are not needed to control any linearization error; they only determine whether the current linearization domain is large enough to contain the certified state and input deviations. Failures such as :inconsistent_at_max_box therefore indicate that the certified ellipsoid requires a larger admissible state or input region than allowed by the current box limits, rather than any issue with the quality of the linear approximation.
# As a result, this toy example is useful for debugging the LMI formulation and controller synthesis independently of nonlinear approximation effects. Any certification failure can be attributed to feasibility, numerical conditioning, input constraints, or box-consistency checks, but not to linearization error.
ellip_opts = EB.EllipsoidalBackwardOptions(;
    maxδx = 1e6, # Upper bound on predecessor ellipsoid size. Larger makes the LMI easier. 
    maxδu = 1e6, # Upper bound on controller/input deviation (deviation from nominal input). Larger makes the LMI easier.
    λ = 0.3, # Objective tradeoff. Small λ favors larger ellipsoids; large λ favors lower transition cost.
    terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.5^2, # shape matrix Q of the terminal ellipsoid centered at the last trajectory point
    transition_cost = trans_cost,
    linearization_δx = [1.1, 1.0], # State box radius used to compute local affine approximation and Lipschitz bounds. used if not using adaptive boxes
    linearization_δu = [0.5, 0.5], # Input box radius used to compute local affine approximation and Lipschitz bounds. used if not using adaptive boxes
    adaptive_boxes = adaptive_opts,
    use_log_det = false,
)

# So ell is allowed to move away from nominal u, but only within the allowed δu budget.

sdp_optimizer = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

certifier = EB.TrajectoryCertifier(provider, sdp_optimizer, ellip_opts)

# ------------------------------------------------------------
# 6) Modular trajectory-generation + certification optimizer
# ------------------------------------------------------------

tc_optimizer = AB.TrajectoryCertificationOptimizer.Optimizer(combo_gen, certifier)

MOI.set(tc_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

MOI.optimize!(tc_optimizer)

composite_traj = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("trajectory"))
controller = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("controller"))
success = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("success"))

println("Trajectory generation + certification success: ", success)
println(
    "Trajectory generation + certification time: ",
    MOI.get(tc_optimizer, MOI.SolveTimeSec()),
)

cert_result = EB.get_result(certifier)

ellipsoids = cert_result.lmi_data.ellipsoids

@show length(composite_traj.x.seq)
@show length(ellipsoids)
@show cert_result.success
@show cert_result.failed_k

# ------------------------------------------------------------
# 7) Plot
# ------------------------------------------------------------

function plot_ellipsoid_chain!(
    fig,
    cert_result;
    max_ellipsoids = 100,
    label_prefix = "Backward ellipsoid",
)
    cert_result === nothing && return fig
    cert_result.lmi_data === nothing && return fig

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

function input_image_ellipsoid(E, κ)
    K = Matrix{Float64}(κ.A)
    b = collect(Float64, κ.c)

    c_u = K * collect(Float64, LazySets.center(E)) + b
    Q_u = K * Matrix{Float64}(LazySets.shape_matrix(E)) * K'

    return LazySets.Ellipsoid(c_u, UT._symmetrize(Q_u))
end

function sampled_indices(n::Int, max_items::Int)
    n <= 0 && return Int[]
    max_items <= 0 && return Int[]

    m = min(max_items, n)
    m == 1 && return [1]

    return unique(round.(Int, range(1, n; length = m)))
end

function plot_input_ellipsoids!(fig, cert_result, Uset; max_ellipsoids = 20)
    cert_result === nothing && return fig

    steps = cert_result.steps
    isempty(steps) && return fig

    valid_steps = filter(s -> s.ellipsoid !== nothing && s.kappa !== nothing, steps)
    isempty(valid_steps) && return fig

    idxs = sampled_indices(length(valid_steps), max_ellipsoids)

    plot!(fig, Uset; alpha = 0.2, label = "Input set")

    for (j, idx) in enumerate(idxs)
        step = valid_steps[idx]

        Uell = input_image_ellipsoid(step.ellipsoid, step.kappa)

        plot!(
            fig,
            Uell;
            linewidth = 1.5,
            linestyle = :dash,
            alpha = 0.6,
            label = j == 1 ? "κ(Eₖ)" : "",
        )
    end

    return fig
end

fig = plot(; aspect_ratio = :equal)

plot!(fig, concrete_problem; aspect_ratio = :equal)

plot_ellipsoid_chain!(fig, cert_result; max_ellipsoids = 100)

if composite_traj !== nothing
    plot!(
        fig,
        composite_traj.x;
        color = :green,
        dims = [1, 2],
        label = "Generated trajectory",
    )
end

display(fig)

# fig_u = plot(; aspect_ratio = :equal)
# # The controller is not forced to interpolate the nominal input.
# # The nomainal input could be outside the ellipdsoidal input

# # If the nominal input lies near the boundary of U, the certification
# # problem may become infeasible because the controller image κ(Eₖ) must
# # remain entirely inside U. The optimizer may therefore shift the
# # certified controller away from the nominal input to recover feasibility.

# max_ellipsoids = 100

# valid_steps = filter(s -> s.ellipsoid !== nothing && s.kappa !== nothing, cert_result.steps)
# idxs = sampled_indices(length(valid_steps), max_ellipsoids)

# plot_input_ellipsoids!(
#     fig_u,
#     cert_result,
#     _U_extended;
#     max_ellipsoids = max_ellipsoids,
# )

# plot!(
#     fig_u,
#     _U_;
#     alpha = 0.2,
#     label = "Input set",
# )

# scatter!(
#     fig_u,
#     [composite_traj.u.seq[valid_steps[i].k][1] for i in idxs],
#     [composite_traj.u.seq[valid_steps[i].k][2] for i in idxs];
#     label = "Nominal inputs",
#     markersize = 3,
# )

# display(fig_u)

for s in cert_result.steps
    println("k=", s.k, " status=", s.status)
end

for s in cert_result.steps
    if s.ellipsoid !== nothing
        eigs = LA.eigvals(Matrix(UT.get_quadratic_form(s.ellipsoid)))
        println("k=", s.k, " min eig=", minimum(eigs), " max eig=", maximum(eigs))
    end
end
