#!/usr/bin/env julia

"""
Per-stage baseline timings for the trajectory-generation + certification pipeline
(plan.md §7, P0). Each stage of the pipeline is timed separately on deterministic
configurations so later phases can prove their speedups against this table:

  seed      — uniform-grid abstraction build + seed trajectory rollout
  mppi      — MPPI refinement of the seed (reports rollout throughput)
  provider  — repeated `build_affine_approximation` calls (cold vs warm exposes the
              per-call symbolic differentiation + compilation cost)
  chain     — ellipsoidal backward certification of the MPPI trajectory (per-step SDP
              cost, Clarabel)

Usage:
    julia --project=bench bench/trajectory_pipeline.jl [filter]

    filter : only run systems whose name contains this substring
             (e.g. `integrator`); default runs integrator + pendulum.

Output:
    bench/results/trajectory_pipeline.csv
"""

using Printf
using Statistics
using Random
using StaticArrays
import MathOptInterface as MOI
import LazySets
import LinearAlgebra as LA
using Symbolics
import Clarabel
using JuMP: optimizer_with_attributes
import CSV
import DataFrames

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const AB = DI.Optim.Abstraction
const EB = AB.EllipsoidalTrajectoryCertifier

const PROBLEMS = joinpath(dirname(dirname(pathof(Dionysos))), "problems")
include(joinpath(PROBLEMS, "Integrator", "integrator.jl"))
include(joinpath(PROBLEMS, "Pendulum", "simple_pendulum.jl"))

# One row per (system, stage): wall time plus one stage-specific metric.
function _row(system, stage, time_s; metric = "", value = NaN, status = "ok")
    return (;
        system = system,
        stage = stage,
        time_s = time_s,
        metric = metric,
        value = value,
        status = status,
    )
end

# ------------------------------------------------------------
# Integrator (2-D, ẋ = u): the fast regression case.
# ------------------------------------------------------------

function integrator_rows()
    rows = NamedTuple[]
    Δt = 0.3
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(4.0, 4.0))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
    concrete_system = Integrator.system(; _X_ = _X_, _U_ = _U_)

    _I_ = LazySets.Hyperrectangle(; low = SVector(-1.7, -1.7), high = SVector(-1.6, -1.6))
    target_set = UT.set_union([
        LazySets.Hyperrectangle(; low = SVector(-1.0, 3.0), high = SVector(-0.3, 3.7)),
        LazySets.Hyperrectangle(; low = SVector(1.0, 2.0), high = SVector(3.0, 3.7)),
    ])
    concrete_problem =
        PR.OptimalControlProblem(concrete_system, _I_, target_set, nothing, nothing, 20)

    # --- seed ---
    seed_stats = @timed begin
        opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("concrete_problem"),
            PR.AlternatingSimulationProblem(concrete_system, concrete_system.X),
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("state_grid"),
            MP.GridFree(SVector(-2.0, -2.0), SVector(0.2, 0.2)),
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("input_grid"),
            MP.GridFree(SVector(-1.0, -1.0), SVector(0.5, 0.5)),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), Δt)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("approx_mode"),
            AB.UniformGridAbstraction.GROWTH,
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("jacobian_bound"),
            Integrator.jacobian_bound(),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("n_samples"), 1)
        MOI.set(opt, MOI.Silent(), true)
        MOI.optimize!(opt)

        gen = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
            opt;
            initial_state = SVector(-1.65, -1.65),
            concrete = false,
            nstep = 20,
        )
        AB.set_problem!(gen, concrete_problem)
        AB.generate!(gen)
        AB.get_trajectory(gen)
    end
    seed_traj = seed_stats.value
    push!(rows, _row("integrator", "seed", seed_stats.time))
    seed_traj === nothing && return rows

    # --- mppi ---
    nstep, nsamples, niter = 20, 200, 5
    discrete_problem = PR.discretize_problem(concrete_problem, Δt)
    mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
        rng = Random.MersenneTwister(0),
        seed_trajectory = seed_traj,
        nstep = nstep,
        nsamples = nsamples,
        niter = niter,
        λ = 1.0,
        noise_sampler = (rng, u, k) -> SVector(0.3 * randn(rng), 0.3 * randn(rng)),
        project_input = u -> SVector(clamp(u[1], -1.0, 1.0), clamp(u[2], -1.0, 1.0)),
        trajectory_cost = (problem, tr) -> sum(LA.norm(u)^2 for u in ST.inputs(tr)),
        hard_constraint = false,
    )
    AB.set_problem!(mppi, discrete_problem)
    mppi_stats = @timed AB.generate!(mppi)
    traj = AB.get_trajectory(mppi)
    rollouts = nsamples * niter
    push!(
        rows,
        _row(
            "integrator",
            "mppi",
            mppi_stats.time;
            metric = "rollouts_per_s",
            value = rollouts / mppi_stats.time,
        ),
    )
    traj === nothing && return rows

    # --- provider ---
    Symbolics.@variables x[1:2] u[1:2] w[1:2]
    fsymbolic = [x[1] + Δt * (u[1] + w[1]), x[2] + Δt * (u[2] + w[2])]
    Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))
    provider = ST.SymbolicAffineApproximationProvider(
        fsymbolic,
        collect(x),
        collect(u),
        collect(w),
        [0.0, 0.0],
        ST.format_input_set(_U_),
        ST.format_noise_set(Wset),
    )
    xk = collect(first(ST.states(traj)))
    uk = collect(first(ST.inputs(traj)))
    δx, δu = [0.2, 0.2], [0.5, 0.5]

    cold = @timed ST.build_affine_approximation(provider, xk, uk; δx = δx, δu = δu)
    warm_times = [
        (@timed ST.build_affine_approximation(provider, xk, uk; δx = δx, δu = δu)).time
        for _ in 1:10
    ]
    push!(
        rows,
        _row(
            "integrator",
            "provider_cold",
            cold.time;
            metric = "ms_per_call",
            value = 1e3 * cold.time,
        ),
    )
    push!(
        rows,
        _row(
            "integrator",
            "provider_warm",
            sum(warm_times);
            metric = "ms_per_call",
            value = 1e3 * median(warm_times),
        ),
    )

    # --- chain ---
    adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
        false,
        [0.2, 0.2],
        [0.01, 0.01],
        [1.0, 1.0],
        [0.5, 0.5],
        [0.01, 0.01],
        [1.0, 1.0],
        2.0,
        1.05,
        1,
        1e-8,
        false,
        [1.0],
        :first_consistent,
        true,
    )
    ellip_opts = EB.ChainOptions(;
        maxδx = 30,
        maxδu = 1.0,
        λ = 0.05,
        terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.5^2,
        linearization_δx = [0.2, 0.2],
        linearization_δu = [0.5, 0.5],
        adaptive_boxes = adaptive_opts,
        objective = :trace,
    )
    sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    cert = EB.BackwardCertifier(provider, sdp, ellip_opts)
    AB.set_problem!(cert, concrete_problem)
    AB.set_trajectory!(cert, traj)
    chain_stats = @timed AB.certify!(cert)
    res = EB.get_result(cert)
    nsteps = length(ST.inputs(traj))
    push!(
        rows,
        _row(
            "integrator",
            "chain",
            chain_stats.time;
            metric = "ms_per_step",
            value = 1e3 * chain_stats.time / max(nsteps, 1),
            status = res === nothing ? "no_result" :
                     (res.success ? "certified" : "failed_k=$(res.failed_k)"),
        ),
    )
    return rows
end

# ------------------------------------------------------------
# Simple pendulum swing-up: the realistic nonlinear case
# (configuration mirrors research/TrajectoryCertificationOptimizer/simple_pendulum.jl,
# with a fixed rng so the baseline is reproducible).
# ------------------------------------------------------------

function pendulum_rows()
    rows = NamedTuple[]
    params = SimplePendulum.Params(; l = 1.0, g = 9.81)
    concrete_problem = SimplePendulum.optimal_control_problem(;
        params = params,
        objective = "reachability_up_medium_power",
    )

    periodic_dims = SVector(1)
    periods = SVector(2π)
    periodic_start = SVector(-π)
    wrap = UT.get_periodic_wrapper(periodic_dims, periods; start = periodic_start)

    Δt = 0.1

    # --- seed ---
    seed_stats = @timed begin
        opt = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
        MOI.set(opt, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
        MOI.set(opt, MOI.RawOptimizerAttribute("h"), SVector(3.0 * π / 180.0, 0.05))
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("input_grid"),
            MP.GridFree(SVector(0.0), SVector(0.3)),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("time_step"), Δt)
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("approx_mode"),
            AB.UniformGridAbstraction.GROWTH,
        )
        MOI.set(
            opt,
            MOI.RawOptimizerAttribute("jacobian_bound"),
            SimplePendulum.jacobian_bound(params),
        )
        MOI.set(opt, MOI.RawOptimizerAttribute("use_periodic_mapping"), true)
        MOI.set(opt, MOI.RawOptimizerAttribute("periodic_dims"), periodic_dims)
        MOI.set(opt, MOI.RawOptimizerAttribute("periodic_periods"), periods)
        MOI.set(opt, MOI.RawOptimizerAttribute("periodic_start"), periodic_start)
        MOI.set(opt, MOI.RawOptimizerAttribute("early_stop"), true)
        MOI.set(opt, MOI.Silent(), true)

        gen = AB.OptimizerTrajectoryGenerator.TrajectoryGenerator(
            opt;
            initial_state = SVector(0.0, 0.0),
            concrete = false,
            nstep = 100,
        )
        AB.set_problem!(gen, concrete_problem)
        AB.generate!(gen)
        AB.get_trajectory(gen)
    end
    seed_traj = seed_stats.value
    push!(rows, _row("pendulum", "seed", seed_stats.time))
    seed_traj === nothing && return rows

    # --- mppi ---
    nstep, nsamples, niter = 100, 1000, 20
    discrete_problem = PR.discretize_problem(concrete_problem, Δt; num_substeps = 5)

    distance_to_target(x, S) =
        S isa LazySets.UnionSetArray ?
        minimum(distance_to_target(x, Si) for Si in LazySets.array(S)) :
        LA.norm(x - LazySets.center(S))

    trajectory_cost = function (problem, tr)
        xs = ST.states(tr)
        us = ST.inputs(tr)
        target_set =
            UT.set_in_period(problem.target_set, periodic_dims, periods, periodic_start)
        Xperiod =
            UT.set_in_period(problem.system.X, periodic_dims, periods, periodic_start)
        best = minimum(distance_to_target(wrap(x), target_set) for x in xs)
        hit_idx = findfirst(x -> wrap(x) ∈ target_set, xs)
        bonus = hit_idx === nothing ? 0.0 : -1000.0 / hit_idx
        violation = sum(x -> wrap(x) ∈ Xperiod ? 0.0 : 1000.0, xs)
        control = sum(LA.norm(u)^2 for u in us)
        return 100.0 * best + 0.01 * control + violation + bonus
    end

    mppi = AB.MPPITrajectoryGenerator.TrajectoryGenerator(;
        rng = Random.MersenneTwister(0),
        seed_trajectory = seed_traj,
        nstep = nstep,
        nsamples = nsamples,
        niter = niter,
        λ = 1.0,
        noise_sampler = (rng, u, k) -> SVector(0.5 * randn(rng)),
        project_input = u -> SVector(clamp(u[1], -7.0, 7.0)),
        trajectory_cost = trajectory_cost,
        wrap_state = (problem, x) -> wrap(x),
        hard_constraint = false,
    )
    AB.set_problem!(mppi, discrete_problem)
    mppi_stats = @timed AB.generate!(mppi)
    traj = AB.get_trajectory(mppi)
    iters = AB.MPPITrajectoryGenerator.get_diagnostics(mppi).iterations
    rollouts = nsamples * iters
    push!(
        rows,
        _row(
            "pendulum",
            "mppi",
            mppi_stats.time;
            metric = "rollouts_per_s",
            value = rollouts / mppi_stats.time,
        ),
    )
    traj === nothing && return rows

    # --- provider (symbolic RK4 discretization, as in the research driver) ---
    Symbolics.@variables θ ω τ w1 w2 T
    xsym = [θ, ω]
    usym = [τ]
    wsym = [w1, w2]
    f_cont(xloc, uloc) =
        [xloc[2], -(params.g / params.l) * Symbolics.sin(xloc[1]) + uloc[1]]
    f_disc = ST.runge_kutta4(f_cont, xsym, usym, T, 1)
    fsymbolic = Symbolics.substitute([f_disc[1] + w1, f_disc[2] + w2], Dict(T => Δt))
    Wset = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0))
    _U_cert_ = LazySets.Hyperrectangle(; low = SVector(-10.5), high = SVector(10.5))
    provider = ST.SymbolicAffineApproximationProvider(
        fsymbolic,
        xsym,
        usym,
        wsym,
        [0.0, 0.0],
        ST.format_input_set(_U_cert_),
        ST.format_noise_set(Wset),
    )
    xk = collect(first(ST.states(traj)))
    uk = collect(first(ST.inputs(traj)))
    δx, δu = [0.2, 0.4], [1.0]

    cold = @timed ST.build_affine_approximation(provider, xk, uk; δx = δx, δu = δu)
    warm_times = [
        (@timed ST.build_affine_approximation(provider, xk, uk; δx = δx, δu = δu)).time
        for _ in 1:10
    ]
    push!(
        rows,
        _row(
            "pendulum",
            "provider_cold",
            cold.time;
            metric = "ms_per_call",
            value = 1e3 * cold.time,
        ),
    )
    push!(
        rows,
        _row(
            "pendulum",
            "provider_warm",
            sum(warm_times);
            metric = "ms_per_call",
            value = 1e3 * median(warm_times),
        ),
    )

    # --- chain ---
    adaptive_opts = EB.AdaptiveLinearizationBoxOptions(
        true,
        [0.05, 0.10],
        [0.005, 0.005],
        [1.8, 2.5],
        [0.25],
        [0.01],
        [10.0],
        1.5,
        1.05,
        30,
        1e-8,
        false,
        [0.75, 1.0, 1.25, 1.5, 2.0],
        :max_volume,
        true,
    )
    ellip_opts = EB.ChainOptions(;
        maxδx = 1.5,
        maxδu = 3.0,
        λ = 0.001,
        terminal_shape = Matrix{Float64}(LA.I, 2, 2) * 0.8^2,
        linearization_δx = [0.2, 0.4],
        linearization_δu = [1.0],
        adaptive_boxes = adaptive_opts,
        objective = :trace,
    )
    sdp = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    cert = EB.BackwardCertifier(provider, sdp, ellip_opts)
    AB.set_problem!(cert, concrete_problem)
    AB.set_trajectory!(cert, traj)
    chain_stats = @timed AB.certify!(cert)
    res = EB.get_result(cert)
    nsteps = length(ST.inputs(traj))
    push!(
        rows,
        _row(
            "pendulum",
            "chain",
            chain_stats.time;
            metric = "ms_per_step",
            value = 1e3 * chain_stats.time / max(nsteps, 1),
            status = res === nothing ? "no_result" :
                     (res.success ? "certified" : "failed_k=$(res.failed_k)"),
        ),
    )
    return rows
end

# ------------------------------------------------------------
# Driver
# ------------------------------------------------------------

function print_table(rows)
    println("\n" * repeat("=", 92))
    println("TRAJECTORY PIPELINE BASELINE (plan.md P0)")
    println(repeat("=", 92))
    @printf(
        "%-12s %-14s %10s   %-16s %12s   %s\n",
        "System",
        "Stage",
        "Time (s)",
        "Metric",
        "Value",
        "Status"
    )
    println(repeat("-", 92))
    for r in rows
        @printf(
            "%-12s %-14s %10.3f   %-16s %12.3f   %s\n",
            r.system,
            r.stage,
            r.time_s,
            r.metric,
            r.value,
            r.status
        )
    end
    println(repeat("=", 92))
    return nothing
end

function main()
    filter = isempty(ARGS) ? "" : lowercase(ARGS[1])
    systems = [("integrator", integrator_rows), ("pendulum", pendulum_rows)]

    rows = NamedTuple[]
    for (name, f) in systems
        (isempty(filter) || occursin(filter, name)) || continue
        println("→ $name ...")
        try
            append!(rows, f())
        catch e
            msg = sprint(showerror, e)
            println("  FAILED: ", first(split(msg, '\n')))
            push!(rows, _row(name, "FAILED", NaN; status = replace(msg, '\n' => " | ")))
        end
        total = sum(r.time_s for r in rows if r.system == name && !isnan(r.time_s))
        push!(rows, _row(name, "total", total))
    end

    print_table(rows)

    dir = joinpath(@__DIR__, "results")
    isdir(dir) || mkpath(dir)
    path = joinpath(dir, "trajectory_pipeline.csv")
    CSV.write(
        path,
        DataFrames.DataFrame(;
            system = [r.system for r in rows],
            stage = [r.stage for r in rows],
            time_s = [r.time_s for r in rows],
            metric = [r.metric for r in rows],
            value = [r.value for r in rows],
            status = [r.status for r in rows],
        ),
    )
    println("Wrote $path")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
