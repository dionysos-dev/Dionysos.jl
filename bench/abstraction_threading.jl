#!/usr/bin/env julia

"""
Dionysos abstraction-build threading benchmark
==============================================
Measures the speed-up of the *threaded* grid abstraction backend
(`Symbolic.ThreadedBackend`) over the sequential one (`Symbolic.SequentialBackend`)
when building the transition relation with
`Symbolic.compute_abstract_system_from_concrete_system!`.

Usage (run with N threads in the current process):
    julia -t <nthreads> --project=bench bench/abstraction_threading.jl [trials] [resolutions] [dt] [du]

Usage (orchestrate several thread counts from one launch):
    julia --project=bench bench/abstraction_threading.jl [trials] [resolutions] [dt] [du] [threads]

Arguments:
    trials       : repetitions per method (default: 5)
    resolutions  : comma-separated grid resolutions, each an isotropic grid (default: 15).
                   A single N or the triplet N,N,N are equivalent.
    dt           : discretization time step (default: 0.1)
    du           : uniform input grid step (default: 0.5)
    threads      : optional comma-separated thread counts to orchestrate (e.g. 1,4,8)
"""

import LazySets
using Dionysos
using MathematicalSystems
using Statistics
using Printf
using StaticArrays
using LinearAlgebra

# Keep BLAS single-threaded so the measured speed-up reflects the abstraction
# backend, not nested BLAS parallelism.
try
    LinearAlgebra.BLAS.set_num_threads(1)
catch
end

const MP = Dionysos.Mapping
const UT = Dionysos.Utils
const ST = Dionysos.System
const SY = Dionysos.Symbolic
const MS = MathematicalSystems

parse_list(str::AbstractString, ::Type{T}) where {T} =
    [parse(T, s) for s in split(str, ",") if !isempty(s)]

function parse_args()
    a = ARGS
    trials = length(a) >= 1 ? parse(Int, a[1]) : 5
    resolutions = length(a) >= 2 ? parse_list(a[2], Int) : [15]
    # Collapse the triplet N,N,N into single N (users often type 25,25,25 per dimension).
    if length(resolutions) == 3 && all(==(resolutions[1]), resolutions)
        resolutions = [resolutions[1]]
    end
    dt = length(a) >= 3 ? parse(Float64, a[3]) : 0.1
    du = length(a) >= 4 ? parse(Float64, a[4]) : 0.5
    threads_list = length(a) >= 5 ? parse_list(a[5], Int) : [Threads.nthreads()]
    return trials, resolutions, dt, du, threads_list
end

"""Build (state_mapping, input_mapping, continuous_system) for a simple 3D linear system."""
function build_test_system(; n_per_dim::Int, input_step::Float64)
    # State domain [0,1]^3 on a uniform grid.
    lb = SVector(0.0, 0.0, 0.0)
    ub = SVector(1.0, 1.0, 1.0)
    h = (ub - lb) ./ (n_per_dim - 1)
    Xmap = MP.ExplicitGridMapping(MP.GridFree(lb, h))
    MP.cover!(Xmap, UT.box(lb, ub), MP.OUTER)

    # Input domain [-1,1]^3 on a uniform grid.
    lb_u = SVector(-1.0, -1.0, -1.0)
    ub_u = SVector(1.0, 1.0, 1.0)
    h_u = SVector(input_step, input_step, input_step)
    Umap = MP.ExplicitGridMapping(MP.GridFree(lb_u, h_u))
    MP.cover!(Umap, UT.box(lb_u, ub_u), MP.OUTER)

    # Continuous dynamics dx/dt = A x + B u.
    A = @SMatrix [
        0.0 1.0 0.0
        0.0 0.0 1.0
        -1.0 -1.0 -1.0
    ]
    B = @SMatrix [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    F_sys(x, u) = A * x + B * u
    continuous_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(F_sys, 3, 3, nothing, nothing)
    return Xmap, Umap, continuous_system
end

"""Build three discrete-time approximations: Centered, OverApprox, GrowthBound."""
function build_approximations(continuous_system, tstep)
    cont_center = ST.ContinuousTimeCenteredSimulation(continuous_system)
    centered = ST.discretize(cont_center, tstep)

    discrete_system = ST.discretize_continuous_system(continuous_system, tstep)

    function simple_over_approx(elem, u)
        center = LazySets.center(elem)
        radius = LazySets.radius_hyperrectangle(elem)
        new_radius = radius * 1.1 .+ 0.01
        new_center = MS.mapping(discrete_system)(center, u)
        return UT.box(new_center - new_radius, new_center + new_radius)
    end
    over = ST.DiscreteTimeOverApproximationMap(discrete_system, simple_over_approx)

    growth_bound_map(r, u) = @inbounds (r * (1.0 + 0.1 * LinearAlgebra.norm(u)) .+ 0.01)
    growth = ST.DiscreteTimeGrowthBound(discrete_system, growth_bound_map)

    return Dict(
        "CenteredSimulation" => centered,
        "OverApproximation" => over,
        "GrowthBound" => growth,
    )
end

# Sequential vs. threaded execution backends compared by the benchmark.
const SEQUENTIAL_BACKEND = SY.SequentialBackend()
const THREADED_BACKEND = SY.ThreadedBackend()

function _build!(symmodel_builder, approx_obj, backend)
    return SY.compute_abstract_system_from_concrete_system!(
        symmodel_builder(),
        approx_obj;
        execution_backend = backend,
    )
end

# Benchmark a single approximation method (sequential vs. threaded).
function benchmark_method(
    method_name::String,
    symmodel_builder,
    approx_obj,
    trials::Int;
    warmups::Int = 1,
)
    println("  Benchmarking: $method_name")

    # Warm-up (excluded from metrics) to remove compilation / first-run effects.
    for _ in 1:warmups
        _build!(symmodel_builder, approx_obj, SEQUENTIAL_BACKEND)
        _build!(symmodel_builder, approx_obj, THREADED_BACKEND)
    end

    seq_results = Vector{NTuple{4, Any}}(undef, trials)
    first_seq_pool = missing
    for i in 1:trials
        result = @timed _build!(symmodel_builder, approx_obj, SEQUENTIAL_BACKEND)
        pool = hasproperty(result, :gcstats) ? result.gcstats.poolalloc : missing
        i == 1 && (first_seq_pool = pool)
        seq_results[i] = (result.time, result.bytes, result.gctime, pool)
    end

    mt_results = Vector{NTuple{4, Any}}(undef, trials)
    first_mt_pool = missing
    for i in 1:trials
        result = @timed _build!(symmodel_builder, approx_obj, THREADED_BACKEND)
        pool = hasproperty(result, :gcstats) ? result.gcstats.poolalloc : missing
        i == 1 && (first_mt_pool = pool)
        mt_results[i] = (result.time, result.bytes, result.gctime, pool)
    end

    seq_times = getindex.(seq_results, 1)
    mt_times = getindex.(mt_results, 1)
    seq_allocs = getindex.(seq_results, 2)
    mt_allocs = getindex.(mt_results, 2)
    seq_pool = [r[4] for r in seq_results if r[4] !== missing]
    mt_pool = [r[4] for r in mt_results if r[4] !== missing]

    seq_mean_time = mean(seq_times)
    mt_mean_time = mean(mt_times)
    seq_best_time = minimum(seq_times)
    mt_best_time = minimum(mt_times)

    return Dict(
        "method" => method_name,
        "seq_mean_time" => seq_mean_time,
        "seq_std_time" => std(seq_times),
        "mt_mean_time" => mt_mean_time,
        "mt_std_time" => std(mt_times),
        "seq_best_time" => seq_best_time,
        "mt_best_time" => mt_best_time,
        "speedup_mean" => seq_mean_time / mt_mean_time,
        "speedup_best" => seq_best_time / mt_best_time,
        "seq_mean_alloc" => mean(seq_allocs),
        "mt_mean_alloc" => mean(mt_allocs),
        "seq_min_alloc" => minimum(seq_allocs),
        "mt_min_alloc" => minimum(mt_allocs),
        "alloc_ratio_mean" => mean(mt_allocs) / mean(seq_allocs),
        "seq_pool_mean" => isempty(seq_pool) ? missing : mean(seq_pool),
        "mt_pool_mean" => isempty(mt_pool) ? missing : mean(mt_pool),
        "seq_first_pool" => first_seq_pool,
        "mt_first_pool" => first_mt_pool,
        "trials" => trials,
        "warmups" => warmups,
    )
end

function format_bytes(bytes::Real)
    units = ("B", "KB", "MB", "GB", "TB")
    value = float(bytes)
    i = 1
    while value >= 1024 && i < length(units)
        value /= 1024
        i += 1
    end
    return @sprintf("%.2f %s", value, units[i])
end

function print_results_table(
    all_results,
    nthreads::Int,
    n_per_dim::Int,
    dt::Float64,
    du::Float64,
)
    isempty(all_results) && (println("(no results)"); return nothing)
    println("\n" * repeat("=", 108))
    println("ABSTRACTION-BUILD THREADING BENCHMARK")
    println(repeat("=", 108))
    println(
        @sprintf(
            "Threads: %d | Grid: %d^3 (%d states) | dt: %.4f | du: %.4f | trials: %d (warmups: %d)",
            nthreads,
            n_per_dim,
            n_per_dim^3,
            dt,
            du,
            all_results[1]["trials"],
            all_results[1]["warmups"]
        )
    )
    println(repeat("-", 108))
    @printf(
        "%-18s %11s %9s %11s %9s %9s %13s %13s %12s\n",
        "Method",
        "μ_seq (s)",
        "σ_seq",
        "μ_mt (s)",
        "σ_mt",
        "Speedup",
        "Alloc Seq",
        "Alloc MT",
        "MT/Seq"
    )
    println(repeat("-", 108))
    for r in all_results
        @printf(
            "%-18s %11.4f %9.4f %11.4f %9.4f %9.2f %13s %13s %12.2f\n",
            r["method"],
            r["seq_mean_time"],
            r["seq_std_time"],
            r["mt_mean_time"],
            r["mt_std_time"],
            r["speedup_mean"],
            format_bytes(r["seq_mean_alloc"]),
            format_bytes(r["mt_mean_alloc"]),
            r["alloc_ratio_mean"]
        )
    end
    println(repeat("=", 108))
    return nothing
end

function run_single_config(trials, n_per_dim, dt, du)
    nthreads = Threads.nthreads()
    Xmap, Umap, csys = build_test_system(; n_per_dim = n_per_dim, input_step = du)
    approximations = build_approximations(csys, dt)
    symmodel_builder = () -> SY.SymbolicModelList(Xmap, Umap)
    all_results = Dict{String, Any}[]
    for (method_name, approx_obj) in sort(collect(approximations); by = x -> x[1])
        try
            push!(
                all_results,
                merge(
                    benchmark_method(method_name, symmodel_builder, approx_obj, trials),
                    Dict("n_per_dim" => n_per_dim, "threads" => nthreads),
                ),
            )
        catch e
            @warn "Benchmark failed" method_name error = e
        end
    end
    print_results_table(all_results, nthreads, n_per_dim, dt, du)
    return all_results
end

function orchestrate(trials, resolutions, dt, du, threads_list)
    current_threads = Threads.nthreads()
    if length(threads_list) == 1 && threads_list[1] == current_threads
        for ncell in resolutions
            run_single_config(trials, ncell, dt, du)
        end
        return nothing
    end
    script = abspath(@__FILE__)
    for t in threads_list
        if t == current_threads
            println("\n[orchestrate] Using existing process with $t threads")
            for ncell in resolutions
                run_single_config(trials, ncell, dt, du)
            end
        else
            println("\n[orchestrate] Spawning process with $t threads")
            run(
                `julia -t $t --project=bench $script $trials $(join(resolutions, ",")) $dt $du`,
            )
        end
    end
    return nothing
end

run_comprehensive_benchmark() = orchestrate(parse_args()...)

if abspath(PROGRAM_FILE) == @__FILE__
    ts = time()
    try
        run_comprehensive_benchmark()
    catch e
        println("❌ Benchmark failed: $e")
        showerror(stdout, e, catch_backtrace())
        println()
    end
    println(@sprintf("Elapsed: %.2f s", time() - ts))
end
