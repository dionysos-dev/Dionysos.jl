using StaticArrays
using MathematicalSystems
using Dionysos
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

include("../../problems/pendulum/simple_pendulum.jl")

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

l = 1.0
g = 9.81
Δt = 0.1

concrete_problem = SimplePendulum.optimal_control_problem(;
    l = l,
    g = g,
    objective = "reachability_up_medium_power",
)

concrete_system = concrete_problem.system
jacobian_bound = SimplePendulum.jacobian_bound(l, g)

_X_ = concrete_system.X
_U_ = concrete_system.U

# Periodic angular coordinate θ ∈ [-π, π]
periodic_dims = SVector(1)
periods = SVector(2π)
periodic_start = SVector(-π)

wrap = ST.get_periodic_wrapper(periodic_dims, periods; start = periodic_start)

wrap_state = (problem, x) -> wrap(x)

# ------------------------------------------------------------
# 2) Trajectory generator block 1: abstraction-based generator
# ------------------------------------------------------------

hx = SVector(3.0 * π / 180.0, 0.05)

input_grid = MP.GridFree(SVector(0.0), SVector(0.3))

abstraction_optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("concrete_problem"),
    concrete_problem,
)

MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("h"), hx)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("time_step"), Δt)

MOI.set(
    abstraction_optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)

MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("jacobian_bound"), jacobian_bound)

MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("periodic_periods"), periods)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)

MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("early_stop"), true)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("efficient"), true)
MOI.set(abstraction_optimizer, MOI.RawOptimizerAttribute("print_level"), 2)
MOI.set(abstraction_optimizer, MOI.Silent(), true)

trajectory_generator = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
    abstraction_optimizer;
    initial_state = SVector(0.0, 0.0),
    concrete = false,
    nstep = 100,
)

# ------------------------------------------------------------
# 3) Trajectory generator block 2: MPPI generator
# ------------------------------------------------------------

_U_cert_ = UT.HyperRectangle(SVector(-7.0), SVector(7.0))

noise_sampler = function (rng, u, k)
    σ = 0.5
    return SVector(σ * randn(rng))
end

project_input = function (u)
    return SVector(clamp(u[1], -7.0, 7.0))
end

function distance_to_target(x, target_set)
    if target_set isa UT.LazySetUnion
        return minimum(distance_to_target(x, S) for S in target_set.sets)
    else
        return LA.norm(x - UT.get_center(target_set))
    end
end

trajectory_cost = function (problem, traj)
    xs = traj.x.seq
    us = traj.u.seq

    target_set =
        UT.set_in_period(problem.target_set, periodic_dims, periods, periodic_start)

    Xperiod = UT.set_in_period(problem.system.X, periodic_dims, periods, periodic_start)

    target_distances = [distance_to_target(wrap(x), target_set) for x in xs]

    best_target_distance = minimum(target_distances)

    hit_idx = findfirst(x -> wrap(x) ∈ target_set, xs)
    target_bonus = hit_idx === nothing ? 0.0 : -1000.0 / hit_idx

    domain_violation_cost = sum(x -> wrap(x) ∈ Xperiod ? 0.0 : 1000.0, xs)

    control_cost = sum(LA.norm(u)^2 for u in us)

    return 100.0 * best_target_distance +
           0.01 * control_cost +
           domain_violation_cost +
           target_bonus
end

mppi_generator = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
    rng = Random.default_rng(),
    seed_trajectory = nothing,
    nstep = 100,
    nsamples = 1000,
    niter = 20,
    λ = 1.0,
    noise_sampler = noise_sampler,
    project_input = project_input,
    trajectory_cost = trajectory_cost,
    wrap_state = wrap_state,
    hard_constraint = false,
)

# ------------------------------------------------------------
# 4) Composite trajectory generator
# ------------------------------------------------------------

combo_gen = AB.CompositeTrajectoryGenerator.TrajectoryGenerator(
    trajectory_generator,
    mppi_generator;
    tstep = Δt,
    num_substeps = 5,
)

# ------------------------------------------------------------
# 5) Ellipsoidal backward certifier
# ------------------------------------------------------------

const EB = AB.EllipsoidalBackwardTrajectoryCertifier

Symbolics.@variables x[1:2]
Symbolics.@variables u[1:1]
Symbolics.@variables w[1:2]

fsymbolic = [x[1] + Δt * (x[2] + w[1]), x[2] + Δt * (-(g / l) * sin(x[1]) + u[1] + w[2])]

Wformat = UT.HyperRectangle(SVector(0.0, 0.0), SVector(0.0, 0.0))

provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    collect(x),
    collect(u),
    collect(w),
    [0.0, 0.0],
    UT.format_input_set(_U_cert_),
    UT.format_noise_set(Wformat),
)

adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    false,
    [0.2, 0.3],
    [0.01, 0.01],
    [1.0, 1.5],
    [0.5],
    [0.01],
    [3.0],
    2.0,
    1.05,
    3,
    1e-8,
    false,
    [1.0],
    :first_consistent,
    true,
)

ellip_opts = EB.EllipsoidalBackwardOptions(;
    maxδx = 2.0,
    maxδu = 3.0,
    λ = 0.05,
    terminal_shape = Matrix{Float64}(LA.I, 2, 2) / 0.3^2,
    linearization_δx = [0.5, 0.8],
    linearization_δu = [2.0],
    adaptive_boxes = adaptive_opts,
)

certifier = EB.TrajectoryCertifier(provider, Clarabel.Optimizer, ellip_opts)

# ------------------------------------------------------------
# 6) Modular trajectory-generation + certification optimizer
# ------------------------------------------------------------

tc_optimizer = AB.TrajectoryCertificationOptimizer.Optimizer(combo_gen, certifier)

MOI.set(tc_optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

MOI.optimize!(tc_optimizer)

pendulum_traj = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("trajectory"))
controller = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("controller"))
success = MOI.get(tc_optimizer, MOI.RawOptimizerAttribute("success"))
cert_result = EB.get_result(certifier)

println("Pendulum trajectory + certification success: ", success)
println("Solve time: ", MOI.get(tc_optimizer, MOI.SolveTimeSec()))

if pendulum_traj !== nothing
    @show length(pendulum_traj.x.seq)
end

if cert_result !== nothing
    @show length(cert_result.lmi_data.ellipsoids)
    @show cert_result.success
    @show cert_result.failed_k
end

# ------------------------------------------------------------
# 7) Plot helpers
# ------------------------------------------------------------

function sampled_indices(n::Int, max_items::Int)
    n <= 0 && return Int[]
    max_items <= 0 && return Int[]

    m = min(max_items, n)
    m == 1 && return [1]

    return unique(round.(Int, range(1, n; length = m)))
end

function away_from_periodic_seam(E; margin = 0.25)
    θ = E.c[1]
    return θ > -π + margin && θ < π - margin
end

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

    idxs = sampled_indices(length(ellipsoids), max_ellipsoids)

    first_label = true

    for idx in idxs
        E = ellipsoids[idx]

        # The current ellipsoidal certifier works in Euclidean coordinates.
        # Do not plot ellipsoids crossing the angular seam θ = ±π.
        away_from_periodic_seam(E) || continue

        plot!(
            fig,
            E;
            label = first_label ? label_prefix : "",
            linewidth = 1.5,
            linestyle = :dash,
            alpha = 0.6,
        )

        first_label = false
    end

    return fig
end

# ------------------------------------------------------------
# 8) Plot state-space trajectory and certificate
# ------------------------------------------------------------

fig = plot(; aspect_ratio = :equal)

Xplot = UT.set_in_period(concrete_problem.system.X, periodic_dims, periods, periodic_start)

Iplot =
    UT.set_in_period(concrete_problem.initial_set, periodic_dims, periods, periodic_start)

Tplot =
    UT.set_in_period(concrete_problem.target_set, periodic_dims, periods, periodic_start)

plot!(fig, Xplot; color = :grey, opacity = 0.15, label = "X")

plot!(fig, Iplot; color = :green, opacity = 0.25, label = "Initial set")

plot!(fig, Tplot; color = :red, opacity = 0.35, label = "Target set")

plot_ellipsoid_chain!(fig, cert_result; max_ellipsoids = 50)

if pendulum_traj !== nothing
    wrapped_xs = [wrap(x) for x in pendulum_traj.x.seq]
    wrapped_traj = ST.Trajectory(wrapped_xs)

    plot!(
        fig,
        wrapped_traj;
        color = :green,
        dims = [1, 2],
        linewidth = 2,
        label = "Generated trajectory",
    )
end

planner_traj = trajectory_generator.trajectory.x
plot!(
    fig,
    planner_traj;
    color = :red,
    dims = [1, 2],
    linewidth = 2,
    label = "Planner trajectory",
)

xlabel!(fig, "θ")
ylabel!(fig, "ω")

display(fig)
