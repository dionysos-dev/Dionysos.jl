#!/usr/bin/env julia

"""
Dionysos Multithreading Benchmark
=================================
Usage (single thread configuration in current process):
    julia -t <nthreads> --project=. bench/multithreading_benchmark.jl \
            [trials] [resolutions] [dt] [du]

Usage (orchestrated across multiple threads):
    julia --project=. bench/multithreading_benchmark.jl [trials] [resolutions] [dt] [du] [threads]

Arguments:
    trials       : number of repetitions per method (default: 5)
    resolutions  : comma-separated list of grid resolutions (each interpreted as an isotropic grid).
                   Passing a single integer N or the triplet N,N,N are equivalent (default: 15)
    dt           : discretization time step (default: 0.1)
    du           : uniform input grid step (default: 0.5)
    threads      : optional comma-separated list of thread counts to orchestrate (e.g. 1,4,8)
"""

using Dionysos
using MathematicalSystems
using Statistics
using Printf
using StaticArrays
using LinearAlgebra

try
    LinearAlgebra.BLAS.set_num_threads(1)
catch
end

const DO = Dionysos.Domain
const UT = Dionysos.Utils
const ST = Dionysos.System

parse_list(str::AbstractString, ::Type{T}) where {T} =
    [parse(T, s) for s in split(str, ",") if !isempty(s)]

function parse_args()
    a = ARGS
    trials = length(a) >= 1 ? parse(Int, a[1]) : 5
    resolutions = length(a) >= 2 ? parse_list(a[2], Int) : [15]
    # Collapse triplet N,N,N into single N (user often types 25,25,25 thinking per dimension)
    if length(resolutions) == 3 && all(==(resolutions[1]), resolutions)
        resolutions = [resolutions[1]]
    end
    dt = length(a) >= 3 ? parse(Float64, a[3]) : 0.1
    du = length(a) >= 4 ? parse(Float64, a[4]) : 0.5
    threads_list = length(a) >= 5 ? parse_list(a[5], Int) : [Threads.nthreads()]
    return trials, resolutions, dt, du, threads_list
end

const SY = Dionysos.Symbolic
const MS = MathematicalSystems

"""Build (state_domain, input_domain, continuous_system) for a simple 3D linear system."""
function build_test_system(; n_per_dim::Int, input_step::Float64)
    # State domain [0,1]^3 (uniform grid)
    lb = SVector(0.0, 0.0, 0.0)
    ub = SVector(1.0, 1.0, 1.0)
    h = (ub - lb) ./ (n_per_dim - 1)
    Xgrid = DO.GridFree(lb, h)
    Xfull = DO.DomainList(Xgrid)
    DO.add_set!(Xfull, UT.HyperRectangle(lb, ub), DO.OUTER)

    # Input domain [-1,1]^3 (uniform step)
    lb_u = SVector(-1.0, -1.0, -1.0)
    ub_u = SVector(1.0, 1.0, 1.0)
    h_u = SVector(input_step, input_step, input_step)
    Ugrid = DO.GridFree(lb_u, h_u)
    Ufull = DO.DomainList(Ugrid)
    DO.add_set!(Ufull, UT.HyperRectangle(lb_u, ub_u), DO.OUTER)

    # Continuous dynamics dx/dt = A x + B u
    A = @SMatrix [
        0.0 1.0 0.0;
        0.0 0.0 1.0;
        -1.0 -1.0 -1.0
    ]
    B = @SMatrix [
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0
    ]
    F_sys(x, u) = A * x + B * u
    continuous_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(F_sys, 3, 3, nothing, nothing)
    return Xfull, Ufull, continuous_system
end

"""Build three discrete-time approximations: Centered, OverApprox, GrowthBound."""
function build_approximations(continuous_system, tstep)
    cont_center = ST.ContinuousTimeCenteredSimulation(continuous_system)
    centered = ST.discretize(cont_center, tstep)

    discrete_system = ST.discretize_continuous_system(continuous_system, tstep)

    function simple_over_approx(elem, u)
        center = UT.get_center(elem)
        radius = UT.get_r(elem)
        new_radius = radius * 1.1 .+ 0.01
        new_center = MS.mapping(discrete_system)(center, u)
        return UT.HyperRectangle(new_center - new_radius, new_center + new_radius)
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

# Benchmark a single method
function benchmark_method(
    method_name::String,
    symmodel_builder,
    approx_obj,
    trials::Int;
    warmups::Int = 1,
)
    println("  Benchmarking: $method_name")

    # Warm-up (excluded from metrics) to remove compilation / first-run effects
    for _ in 1:warmups
        SY.compute_abstract_system_from_concrete_system!(
            symmodel_builder(),
            approx_obj;
            verbose = false,
            threaded = false,
        )
        SY.compute_abstract_system_from_concrete_system!(
            symmodel_builder(),
            approx_obj;
            verbose = false,
            threaded = true,
        )
    end

    seq_results = Vector{NTuple{4, Any}}(undef, trials)
    first_seq_pool = missing
    for i in 1:trials
        symmodel = symmodel_builder()
        result = @timed SY.compute_abstract_system_from_concrete_system!(
            symmodel,
            approx_obj;
            verbose = false,
            threaded = false,
        )
        pool = hasproperty(result, :gcstats) ? result.gcstats.poolalloc : missing
        if i == 1
            ;
            first_seq_pool = pool;
        end
        seq_results[i] = (result.time, result.bytes, result.gctime, pool)
    end

    mt_results = Vector{NTuple{4, Any}}(undef, trials)
    first_mt_pool = missing
    for i in 1:trials
        symmodel = symmodel_builder()
        result = @timed SY.compute_abstract_system_from_concrete_system!(
            symmodel,
            approx_obj;
            verbose = false,
            threaded = true,
        )
        pool = hasproperty(result, :gcstats) ? result.gcstats.poolalloc : missing
        if i == 1
            ;
            first_mt_pool = pool;
        end
        mt_results[i] = (result.time, result.bytes, result.gctime, pool)
    end

    seq_times = getindex.(seq_results, 1)
    mt_times = getindex.(mt_results, 1)
    seq_allocs = getindex.(seq_results, 2)
    mt_allocs = getindex.(mt_results, 2)
    seq_pool = [r[4] for r in seq_results if r[4] !== missing]
    mt_pool = [r[4] for r in mt_results if r[4] !== missing]

    seq_mean_time = mean(seq_times);
    mt_mean_time = mean(mt_times)
    seq_std_time = std(seq_times);
    mt_std_time = std(mt_times)
    seq_best_time = minimum(seq_times);
    mt_best_time = minimum(mt_times)

    seq_mean_alloc = mean(seq_allocs);
    mt_mean_alloc = mean(mt_allocs)
    seq_min_alloc = minimum(seq_allocs);
    mt_min_alloc = minimum(mt_allocs)

    speedup_mean = seq_mean_time / mt_mean_time
    speedup_best = seq_best_time / mt_best_time  # retained (not displayed)
    alloc_ratio_mean = mt_mean_alloc / seq_mean_alloc
    alloc_ratio_best = mt_min_alloc / seq_min_alloc

    seq_pool_mean = isempty(seq_pool) ? missing : mean(seq_pool)
    mt_pool_mean = isempty(mt_pool) ? missing : mean(mt_pool)

    return Dict(
        "method" => method_name,
        "seq_mean_time" => seq_mean_time,
        "seq_std_time" => seq_std_time,
        "mt_mean_time" => mt_mean_time,
        "mt_std_time" => mt_std_time,
        "seq_best_time" => seq_best_time,
        "mt_best_time" => mt_best_time,
        "speedup_mean" => speedup_mean,
        "speedup_best" => speedup_best,
        "seq_mean_alloc" => seq_mean_alloc,
        "mt_mean_alloc" => mt_mean_alloc,
        "seq_min_alloc" => seq_min_alloc,
        "mt_min_alloc" => mt_min_alloc,
        "alloc_ratio_mean" => alloc_ratio_mean,
        "alloc_ratio_best" => alloc_ratio_best,
        "seq_pool_mean" => seq_pool_mean,
        "mt_pool_mean" => mt_pool_mean,
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
    println("\n" * repeat("=", 108))
    println("MULTITHREADING BENCHMARK RESULTS")
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
    println(repeat("-", 108))
    println("Pool allocation counts (first measured run per method):")
    for r in all_results
        seqp = r["seq_first_pool"] === missing ? "-" : string(r["seq_first_pool"])
        mtp = r["mt_first_pool"] === missing ? "-" : string(r["mt_first_pool"])
        println(@sprintf("  %-18s seq=%s  mt=%s", r["method"], seqp, mtp))
    end
    println(repeat("=", 108))
    return nothing
end

function run_single_config(trials, n_per_dim, dt, du)
    nthreads = Threads.nthreads()
    Xfull, Ufull, csys = build_test_system(; n_per_dim = n_per_dim, input_step = du)
    approximations = build_approximations(csys, dt)
    symmodel_builder = () -> SY.SymbolicModelList(Xfull, Ufull)
    all_results = Dict{String, Any}[]
    for (method_name, approx_obj) in sort(collect(approximations); by = x->x[1])
        try
            push!(
                all_results,
                merge(
                    benchmark_method(method_name, symmodel_builder, approx_obj, trials),
                    Dict("n_per_dim"=>n_per_dim, "threads"=>nthreads),
                ),
            )
        catch e
            @warn "Benchmark failed" method_name error=e
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
            cmd = `julia -t $t --project=. $script $trials $(join(resolutions,",")) $dt $du`
            run(cmd)
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
        showerror(stdout, e, catch_backtrace());
        println()
    end
    println(@sprintf("Elapsed: %.2f s", time() - ts))
end
