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

include("../../problems/Pendulum/simple_pendulum.jl")

# ------------------------------------------------------------
# 1) Problem
# ------------------------------------------------------------

params = SimplePendulum.Params(; l = 1.0, g = 9.81)

concrete_problem = SimplePendulum.optimal_control_problem(;
    params = params,
    objective = "reachability_up_medium_power",
)
concrete_system = concrete_problem.system
jacobian_bound = SimplePendulum.jacobian_bound(params)

# ------------------------------------------------------------
# 2) Trajectory generator block 1: abstraction-based generator
# ------------------------------------------------------------

_X_ = concrete_system.X
_U_ = concrete_system.U

# Periodic angular coordinate θ ∈ [-π, π]
periodic_dims = SVector(1)
periods = SVector(2π)
periodic_start = SVector(-π)

wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = periodic_start)
wrap_state = (problem, x) -> wrap(x)

Δt = 0.1
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

k=1.5
_U_cert_ = UT.box(SVector(-7.0)*k, SVector(7.0)*k)

noise_sampler = function (rng, u, k)
    σ = 0.5
    return SVector(σ * randn(rng))
end

project_input = function (u)
    return SVector(clamp(u[1], -7.0, 7.0))
end

function distance_to_target(x, target_set)
    if target_set isa UT.LazySets.UnionSetArray
        return minimum(distance_to_target(x, S) for S in UT.LazySets.array(target_set))
    else
        return LA.norm(x - LazySets.center(target_set))
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
    Δt = Δt,
    num_substeps = 5,
)

# ------------------------------------------------------------
# 5) Ellipsoidal backward certifier
# ------------------------------------------------------------

const EB = AB.EllipsoidalBackwardTrajectoryCertifier

# Symbolics.@variables x[1:2]
# Symbolics.@variables u[1:1]
# Symbolics.@variables w[1:2]

# fsymbolic = [
#     x[1] + Δt * (x[2] + w[1]),
#     x[2] + Δt * (-(params.g / params.l) * Symbolics.sin(x[1]) + u[1] + w[2]),
# ]
# Wformat = UT.box(SVector(0.0, 0.0), SVector(0.0, 0.0))

# provider = ST.SymbolicAffineApproximationProvider(
#     fsymbolic,
#     collect(x),
#     collect(u),
#     collect(w),
#     [0.0, 0.0],
#     ST.format_input_set(_U_cert_),
#     ST.format_noise_set(Wformat),
# )

Symbolics.@variables θ ω τ w1 w2 T

x = [θ, ω]
u = [τ]
w = [w1, w2]

f_cont_expr(xloc, uloc) =
    [xloc[2], -(params.g / params.l) * Symbolics.sin(xloc[1]) + uloc[1]]

f_disc = ST.runge_kutta4(f_cont_expr, x, u, T, 1)

fsymbolic = Symbolics.substitute([f_disc[1] + w1, f_disc[2] + w2], Dict(T => Δt))

Wset = UT.box(SVector(0.0, 0.0), SVector(0.0, 0.0))

provider = ST.SymbolicAffineApproximationProvider(
    fsymbolic,
    x,
    u,
    w,
    [0.0, 0.0],
    ST.format_input_set(_U_cert_),
    ST.format_noise_set(Wset),
)

adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
    true,                  # Enable adaptive box search
    [0.05, 0.10],          # Initial state box radius δx
    [0.005, 0.005],        # Minimum allowed δx
    [1.8, 2.5],            # Maximum allowed δx
    [0.25],                # Initial input box radius δu
    [0.01],                # Minimum allowed δu
    [10.0],                 # Maximum allowed δu
    1.5,                   # Grow factor
    1.05,                  # Shrink factor
    30,                    # Maximum adaptive iterations
    1e-8,                  # Numerical tolerance
    false,                 # Do not use explicit candidate scales
    [0.75, 1.0, 1.25, 1.5, 2.0],  # Candidate scale list
    # (ignored when previous flag=false)
    :max_volume,           # Select feasible candidate with largest volume
    true,                 # Verbose/debug output
)

ellip_opts = EB.EllipsoidalBackwardOptions(;
    maxδx = 1.5,
    maxδu = 3.0,
    λ = 0.001,
    terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.8^2, # shape matrix Q: semi-axes² on the diagonal
    linearization_δx = [0.2, 0.4],
    linearization_δu = [1.0],
    adaptive_boxes = adaptive_opts,
    use_log_det = true,
)

sdp_optimizer = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

certifier = EB.TrajectoryCertifier(provider, sdp_optimizer, ellip_opts)

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
    θ = LazySets.center(E)[1]
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
        # away_from_periodic_seam(E) || continue

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

fig_u = plot(; xlabel = "k", ylabel = "u", legend = :best)

valid_steps = filter(s -> s.ellipsoid !== nothing && s.kappa !== nothing, cert_result.steps)

idxs = sampled_indices(length(valid_steps), 100)

# Plot physical/certification input limits
hline!(
    fig_u,
    [LazySets.low(_U_cert_, 1), LazySets.high(_U_cert_, 1)];
    label = "Certified input bounds",
    linestyle = :dash,
)
# hline!(fig_u, [LazySets.low(_U_, 1), LazySets.high(_U_, 1)]; label = "Problem input bounds", linestyle = :dot)

# Nominal inputs
scatter!(
    fig_u,
    [valid_steps[i].k for i in idxs],
    [pendulum_traj.u.seq[valid_steps[i].k][1] for i in idxs];
    label = "Nominal input",
    markersize = 3,
)

# Controller image intervals κ(Eₖ)
for (j, i) in enumerate(idxs)
    step = valid_steps[i]
    Uell = input_image_ellipsoid(step.ellipsoid, step.kappa)

    c = LazySets.center(Uell)[1]
    r = sqrt(LazySets.shape_matrix(Uell)[1, 1])

    plot!(
        fig_u,
        [step.k, step.k],
        [c - r, c + r];
        linewidth = 2,
        label = j == 1 ? "κ(Eₖ)" : "",
    )
end

display(fig_u)

for s in cert_result.steps
    println("k=", s.k, " status=", s.status)
end

for s in cert_result.steps
    if s.ellipsoid !== nothing
        eigs = LA.eigvals(Matrix(UT.get_quadratic_form(s.ellipsoid)))
        println("k=", s.k, " min eig=", minimum(eigs), " max eig=", maximum(eigs))
    end
end

s = cert_result.steps[end]

@show s.k
@show s.status
@show s.summary.adaptive_box_status
@show s.summary.Xbar_radius
@show s.summary.Ubar_radius
@show s.summary.required_X_radius
@show s.summary.required_U_radius

idx_fail = findfirst(s -> s.status == :infeasible, cert_result.steps)

if idx_fail !== nothing
    s_fail = cert_result.steps[idx_fail]

    @show s_fail.k
    @show s_fail.summary.adaptive_box_status
    @show s_fail.summary.Xbar_radius
    @show s_fail.summary.Ubar_radius
    @show s_fail.summary.required_X_radius
    @show s_fail.summary.required_U_radius
    @show s_fail.summary.provider_summary

    k = s_fail.k

    xk = collect(pendulum_traj.x.seq[k])
    xnext = collect(pendulum_traj.x.seq[k + 1])
    uk = collect(pendulum_traj.u.seq[k])

    @show k
    @show xk
    @show xnext
    @show uk
    @show maximum(abs.(uk))
    @show xnext .- xk
    @show LA.norm(xnext .- xk)
else
    println("No failed step: certification succeeded.")
end

system_plot! = SimplePendulum.system_plot!(; params = params)
Dionysos.animate_trajectory_dashboard(
    system_plot!,
    pendulum_traj.x,
    pendulum_traj.u;
    xdims = (1, 2),      # phase plot θ vs ω
    udims = (1,),        # input over time
    Δt = Δt,
    fps = 5,
    filename = "pendulum_dashboard.mp4",
    xlabel_state = "θ [rad]",
    ylabel_state = "ω [rad/s]",
    xlabel_input = "time [s]",
    ylabel_input = "τ [Nm]",
    xlims_state = (-π, π),
    ylims_state = (-8, 8),
)
