# Shared pendulum swing-up machinery for campaigns C2/C3 (plan.md §6b): the
# demo_pendulum.jl generation stack (abstraction seed → CEM refinement with the
# wrap-aware terminal pull → periodic unwrap) as a per-rng trajectory factory,
# plus the z-frame helpers every certifier config needs (normalized provider,
# z-scaled problem boxes, backward/forward option builders, funnel-area metric).
#
# Generation is cached across configs: the campaign runner reseeds
# `MersenneTwister(1000 + s)` per (config, seed), so the first UInt64 draw is a
# stable per-seed key — every config certifies the SAME trajectory, and the
# comparison is paired.

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const AB = DI.Optim.Abstraction
const EB = AB.EllipsoidalTrajectoryCertifier
const MPPI = AB.MPPITrajectoryGenerator

import LazySets
import LinearAlgebra as LA
import MathOptInterface as MOI
import MathematicalSystems as MS
using StaticArrays
using Random
using Statistics
using Symbolics
import Clarabel
using JuMP: optimizer_with_attributes

const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")
include(joinpath(PROBLEMS, "Pendulum", "simple_pendulum.jl"))

const pend_params = SimplePendulum.Params(; l = 1.0, g = 9.81)
const pend_problem = SimplePendulum.optimal_control_problem(;
    params = pend_params,
    objective = "reachability_up_medium_power_no_obstacle",
)
const Δt = 0.1
const periodic_dims = SVector(1)
const periods = SVector(2π)
const pend_wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = SVector(-π))
const u_max = 4.5
const u_plan = 0.6 * u_max          # input reserve: the certificate gets the rest

# The demo's fixed global normalization (85% of the target θ half-width, ω 0.25).
const T_FIXED = [0.85 * 15.0 * π / 180.0, 0.25]

const sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

# ------------------------------------------------------------
# Generation (demo_pendulum.jl stack, verbatim)
# ------------------------------------------------------------

function _abstraction_seed()
    opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
    MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), pend_problem)
    MOI.set(opt, MOI.RawOptimizerAttribute("h"), SVector(3.0 * π / 180.0, 0.1))
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("input_grid"),
        MP.GridFree(SVector(0.0), SVector(0.3)),
    )
    MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), Δt)
    MOI.set(opt, MOI.RawOptimizerAttribute("approx_mode"), AB.UniformGridAbstraction.GROWTH)
    MOI.set(
        opt,
        MOI.RawOptimizerAttribute("jacobian_bound"),
        SimplePendulum.jacobian_bound(pend_params),
    )
    MOI.set(opt, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
    MOI.set(opt, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
    MOI.set(opt, MOI.RawOptimizerAttribute("periodic_periods"), periods)
    MOI.set(opt, MOI.RawOptimizerAttribute("periodic_start"), SVector(-π))
    MOI.set(opt, MOI.RawOptimizerAttribute("early_stop"), true)
    MOI.set(opt, MOI.Silent(), true)

    seed_gen = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
        opt;
        initial_state = SVector(0.0, 0.0),
        concrete = false,
        nstep = 100,
    )
    AB.set_problem!(seed_gen, pend_problem)
    AB.generate!(seed_gen)
    return AB.get_trajectory(seed_gen)
end
const pend_seed = _abstraction_seed()

const pend_base = PR.discretize_problem(pend_problem, Δt; num_substeps = 5)
const T_split =
    UT.set_in_period(pend_problem.target_set, periodic_dims, periods, SVector(-π))
const pend_discrete = PR.OptimalControlProblem(
    pend_base.system,
    pend_base.initial_set,
    T_split,
    pend_base.state_cost,
    pend_base.transition_cost,
    pend_base.time,
    pend_base.safe_set,
)

const pend_cost = AB.CompositeCost(
    AB.ReachObjectiveCost(T_split; wrap = pend_wrap),
    AB.TerminalPullCost(
        [π, 0.0],
        [0.249, 0.95];
        w = 500.0,
        wrap = pend_wrap,
        periods = [2π, nothing],
    ),
    AB.InputEffortCost(0.001),
    AB.InputSmoothnessCost(; w_du = 0.05, w_ddu = 0.01),
    AB.DomainPenaltyCost(pend_problem.system.X; wrap = pend_wrap),
)

# x-frame lift: periodic unwrap, shifted so the endpoint lands in the target's
# θ-range.
function lift_pendulum(traj)
    lifted = ST.unwrap_trajectory(traj, (1,), (2π,))
    θN = collect(ST.states(lifted))[end][1]
    shift = 2π * round((θN - π) / (2π))
    shift == 0.0 && return lifted
    xs = [SVector(x[1] - shift, x[2]) for x in ST.states(lifted)]
    return ST.Trajectory(xs; inputs = collect(ST.inputs(lifted)))
end

function generate_pendulum(rng)
    mppi = MPPI.TrajectoryGenerator(;
        rng = rng,
        seed_trajectory = pend_seed,
        nstep = 90,
        nsamples = 1000,
        niter = 40,
        noise = MPPI.GaussianMPPINoise(SVector(0.8)),
        project_input = u -> SVector(clamp(u[1], -u_plan, u_plan)),
        cost = pend_cost,
        wrap_state = (problem, x) -> pend_wrap(x),
        update_rule = :cem,
        elite_frac = 0.05,
        antithetic = true,
        stop_on_success = false,
    )
    AB.set_problem!(mppi, pend_discrete)
    stats = @timed AB.generate!(mppi)
    traj = AB.get_trajectory(mppi)
    traj === nothing && return nothing
    return (;
        traj = lift_pendulum(traj),
        gen_success = AB.get_success(mppi),
        gen_time = stats.time,
    )
end

# Two cache layers: in-process across configs, and on disk across campaign
# processes (C2 and C3 draw the same per-seed trajectories). The disk cache is
# regenerable scratch — `campaigns/cache/` is gitignored.
const _PEND_CACHE = Dict{UInt64, Any}()
const _PEND_DISK = joinpath(@__DIR__, "cache")

function _save_pendulum(key, gen)
    isdir(_PEND_DISK) || mkpath(_PEND_DISK)
    open(joinpath(_PEND_DISK, "pend_$key.txt"), "w") do io
        println(io, gen.gen_success ? 1 : 0)
        for x in ST.states(gen.traj)
            println(io, join(x, ","))
        end
        println(io, "#")
        for u in ST.inputs(gen.traj)
            println(io, join(u, ","))
        end
    end
    return nothing
end

function _load_pendulum(key)
    path = joinpath(_PEND_DISK, "pend_$key.txt")
    isfile(path) || return nothing
    lines = readlines(path)
    isep = findfirst(==("#"), lines)
    isep === nothing && return nothing
    xs = [SVector{2}(parse.(Float64, split(l, ","))...) for l in lines[2:(isep - 1)]]
    us = [SVector{1}(parse(Float64, l)) for l in lines[(isep + 1):end]]
    return (;
        traj = ST.Trajectory(xs; inputs = us),
        gen_success = lines[1] == "1",
        gen_time = 0.0,
    )
end

function cached_pendulum(rng)
    key = rand(rng, UInt64)
    return get!(_PEND_CACHE, key) do
        loaded = _load_pendulum(key)
        loaded === nothing || return loaded
        gen = generate_pendulum(rng)
        gen === nothing || _save_pendulum(key, gen)
        return gen
    end
end

# ------------------------------------------------------------
# Certification frame (z = x ./ t, plan.md §4.3)
# ------------------------------------------------------------

function _raw_provider()
    Symbolics.@variables θ ω τ w1 w2 T
    f_cont(xloc, uloc) = collect(SimplePendulum.dynamic(pend_params)(xloc, uloc))
    f_disc = ST.runge_kutta4(f_cont, [θ, ω], [τ], T, 1)
    fsymbolic = Symbolics.substitute([f_disc[1] + w1, f_disc[2] + w2], Dict(T => Δt))
    Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))
    return ST.SymbolicAffineApproximationProvider(
        fsymbolic,
        [θ, ω],
        [τ],
        [w1, w2],
        [0.0, 0.0],
        ST.format_input_set(pend_problem.system.U),
        ST.format_noise_set(Wset),
    )
end
const pend_provider = _raw_provider()

const _ZPROVIDERS = Dict{Vector{Float64}, Any}()
function zprovider(t)
    return get!(
        () -> ST.normalized_symbolic_provider(pend_provider, collect(t)),
        _ZPROVIDERS,
        collect(Float64, t),
    )
end

function zproblem(t)
    return PR.OptimalControlProblem(
        pend_base.system,
        ST.normalize_box(pend_problem.initial_set, t),
        ST.normalize_box(pend_problem.target_set, t),
        nothing,
        nothing,
        pend_base.time,
        nothing,
    )
end
ztraj(traj, t) = ST.normalize_trajectory(traj, t)

# Florentin's `trajectory_std` normalization: t from the per-dimension spread of
# the trajectory itself (floored — a flat dimension must not blow the frame up).
function trajectory_std(traj)
    xs = collect(ST.states(traj))
    return [max(Statistics.std([x[i] for x in xs]), 1e-3) for i in 1:2]
end

# The demo's certified configuration, parameterized on the C2 axes. The
# linearization-box search is fixed at :max_volume so the ablation axes are the
# only moving parts. `check_state_domain = false`: the lifted θ leaves the
# wrapped domain by construction and ω has a generous margin (see the demo).
function pend_backward_options(
    t;
    objective = :maximin,
    remainder_model = :vertices,
    box_search = :max_volume,
)
    adaptive = EB.AdaptiveLinearizationBoxOptions(;
        ΔX_initial = [0.05, 0.10] ./ t,
        ΔX_min = [0.005, 0.005] ./ t,
        ΔX_max = [2.5, 3.5] ./ t,
        ΔU_initial = [0.25],
        ΔU_min = [0.01],
        ΔU_max = [4.5],
        search_scales = [0.75, 1.0, 1.5, 2.0],
        objective = box_search,
    )
    return EB.ChainOptions(;
        maxδx = 12.0,
        maxδu = 3.0,
        λ = 0.001,
        terminal_shape = nothing,
        terminal_shrink = 0.95,
        linearization_δx = [0.05, 0.10] ./ t,
        linearization_δu = [1.0],
        adaptive_boxes = adaptive,
        objective = objective,
        remainder_model = remainder_model,
        check_state_domain = false,
    )
end

# Mode-asymmetric λ and an entry-scale q_min floor — see the C3 script for the
# measured razor-edge (:fixed, λ ≪ 1) and drift/collapse (:free, λ big / floor
# loose) mechanisms. The z-frame floor 1.0 sits just under the full entry's
# smallest z-radius (θ: √2·0.784 ≈ 1.11), mirroring the integrator calibration.
function pend_forward_options(;
    target_mode = :fixed,
    remainder_model = :vertices,
    entry_shape = nothing,
)
    return EB.ForwardOptions(;
        target_mode = target_mode,
        entry_shape = entry_shape,
        maxδu = 3.0,
        λ = target_mode === :fixed ? 0.5 : 0.01,
        q_min = 1.0,
        linearization_δu = [1.0],
        linearization_δx_margin = 1.1,
        check_state_domain = false,
        remainder_model = remainder_model,
    )
end

# The certifier's circumscribed entry (√n·halfwidths), scaled: the forward
# entry ladder — how much of the initial set can a forward tube defend?
function pend_entry_shape(t, scale)
    zI = ST.normalize_box(pend_problem.initial_set, t)
    a = sqrt(2.0) .* scale .* collect(LazySets.radius_hyperrectangle(zI))
    return Matrix(LA.Diagonal(a .^ 2))
end

# Funnel areas in the PHYSICAL frame (V = π·√det(D·Q·D)), comparable across
# normalizations and to the PR sweep report.
function funnel_areas(res, t)
    res === nothing && return Float64[]
    D = LA.Diagonal(collect(Float64, t))
    return [
        π * sqrt(max(0.0, LA.det(D * Matrix(LazySets.shape_matrix(E)) * D))) for
        E in res.lmi_data.ellipsoids
    ]
end
