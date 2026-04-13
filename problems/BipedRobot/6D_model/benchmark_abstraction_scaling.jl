# ==============================================================================
# Benchmark distributed abstraction scaling
# ==============================================================================
#
# Measures separately:
#   - master package load time      (only once per Julia session)
#   - master include time           (only once per Julia session)
#   - worker creation time
#   - worker package load time
#   - worker include time
#   - worker warm-up time
#   - concrete problem build time
#   - optimizer build time
#   - abstraction optimize! wall time
#   - abstraction reported construction time
#   - end-to-end benchmark time
#
# Notes:
#   * Repeats in the same Julia session are cold for workers but warm for master.
#   * If you want a fully cold benchmark including master startup/package loading,
#     launch a fresh Julia process for each run.
# ==============================================================================

using Distributed
using Dates
using Printf
using Statistics

# ------------------------------------------------------------------------------
# Benchmark configuration
# ------------------------------------------------------------------------------

WORKER_COUNTS = [0, 1, 2, 4, 8]
NREPEAT = 3

TSTEP = 0.1
DISTRIBUTED_PARTITION_STRATEGY = :contiguous   # or :roundrobin
USE_THREADED_PER_WORKER = false

PRINT_LEVEL = 1
PROGRESS_UPDATE_INTERVAL = Int(1e2)
PROGRESS_DT = 60.0
SIMPLIFY = 3.0  # grid cell size multiplier for faster abstraction construction during benchmark

SAVE_RESULTS = false
RESULTS_FILE = joinpath(@__DIR__, "benchmark_abstraction_scaling_results.txt")

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

ROBOT_PROBLEM_PATH = joinpath(@__DIR__, "robot_problem.jl")
UTILS_PATH = joinpath(@__DIR__, "utils.jl")
ROBOT_URDF = joinpath(@__DIR__, "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf")

# ------------------------------------------------------------------------------
# Master startup timing
# ------------------------------------------------------------------------------

const MASTER_STARTUP = let
    t_master_packages = @elapsed begin
        using Dionysos
        using MathematicalSystems
        using StaticArrays
        using LinearAlgebra
        using JuMP
        using MathOptInterface
    end

    MOI = MathOptInterface

    const DI = Dionysos
    const UT = DI.Utils
    const ST = DI.System
    const MP = DI.Mapping
    const OP = DI.Optim
    const AB = OP.Abstraction

    t_master_includes = @elapsed begin
        include(ROBOT_PROBLEM_PATH)
        include(UTILS_PATH)
        using .RobotProblem
    end

    (t_master_packages = t_master_packages, t_master_includes = t_master_includes)
end

const MOI = MathOptInterface
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const OP = DI.Optim
const AB = OP.Abstraction

# ------------------------------------------------------------------------------
# Worker reset
# ------------------------------------------------------------------------------

function reset_workers!()
    ps = workers()
    isempty(ps) || rmprocs(ps)
    return nothing
end

# ------------------------------------------------------------------------------
# Worker setup split into phases
# ------------------------------------------------------------------------------

function setup_workers!(nworkers::Int)
    reset_workers!()

    t_addprocs = 0.0
    t_worker_packages = 0.0
    t_worker_includes = 0.0

    if nworkers > 0
        t_addprocs = @elapsed addprocs(nworkers; exeflags="--project=$(Base.active_project())")

        t_worker_packages = @elapsed begin
            @everywhere begin
                using Dionysos
                using MathematicalSystems
                using StaticArrays
                using LinearAlgebra
            end
        end

        t_worker_includes = @elapsed begin
            @everywhere include($ROBOT_PROBLEM_PATH)
            @everywhere include($UTILS_PATH)
            @everywhere using .RobotProblem
        end
    end

    return (
        t_addprocs = t_addprocs,
        t_worker_packages = t_worker_packages,
        t_worker_includes = t_worker_includes,
    )
end

# ------------------------------------------------------------------------------
# Warm-up
# ------------------------------------------------------------------------------

function warmup_workers!(; robot_urdf, tstep)
    isempty(workers()) && return 0.0

    t = @elapsed begin
        futures = [
            remotecall(
                RobotProblem.warmup_robot_problem!,
                p;
                robot_urdf = robot_urdf,
                tstep = tstep,
            ) for p in workers()
        ]
        fetch.(futures)
    end

    return t
end

# ------------------------------------------------------------------------------
# Optimizer builder
# ------------------------------------------------------------------------------

function build_optimizer(;
    concrete_problem,
    state_grid,
    input_grid,
    distributed::Bool,
    distributed_nparts::Int,
    distributed_partition_strategy::Symbol,
    threaded::Bool,
    state_filter = nothing,
    state_input_filter = nothing,
    print_level::Int = PRINT_LEVEL,
)
    optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)

    concrete_system = concrete_problem.system

    MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), state_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("use_implicit_mapping"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("mapping_region"), concrete_system.X)

    MOI.set(optimizer, MOI.RawOptimizerAttribute("input_grid"), input_grid)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_filter"), state_filter)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("state_input_filter"), state_input_filter)

    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("approx_mode"),
        AB.UniformGridAbstraction.CENTER_SIMULATION,
    )

    MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed"), distributed)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("distributed_nparts"), distributed_nparts)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("distributed_partition_strategy"),
        distributed_partition_strategy,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("threaded"), threaded)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("efficient"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("print_level"), print_level)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("progress_update_interval"),
        PROGRESS_UPDATE_INTERVAL,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("progress_dt"), PROGRESS_DT)

    return optimizer
end

# ------------------------------------------------------------------------------
# One benchmark run
# ------------------------------------------------------------------------------

function run_one_benchmark(nworkers::Int)
    distributed = nworkers > 0

    t_end_to_end = @elapsed begin
        # 1) Setup workers
        worker_setup = setup_workers!(nworkers)

        # 2) Warm-up workers
        t_warmup = warmup_workers!(; robot_urdf = ROBOT_URDF, tstep = TSTEP)

        # 3) Build concrete problem
        t_problem_build = @elapsed begin
            global concrete_problem =
                RobotProblem.problem(; robot_urdf = ROBOT_URDF, tstep = TSTEP)
            global concrete_system = concrete_problem.system
        end

        n_state = MathematicalSystems.statedim(concrete_system)
        n_input = MathematicalSystems.inputdim(concrete_system)

        # 4) Build grids
        x0 = SVector{n_state, Float64}(zeros(n_state))
        hx =
            SVector{n_state, Float64}([fill(2π / 180, 3)..., fill(0.15, 3)...]) * SIMPLIFY
        state_grid = MP.GridFree(x0, hx)

        u0 = SVector{n_input, Float64}(zeros(n_input))
        hu = SVector{n_input, Float64}(fill(1.0, n_input))
        input_grid = MP.GridFree(u0, hu)

        state_filter = nothing
        state_input_filter = nothing

        # 5) Build optimizer
        t_optimizer_build = @elapsed begin
            global optimizer = build_optimizer(;
                concrete_problem = concrete_problem,
                state_grid = state_grid,
                input_grid = input_grid,
                distributed = distributed,
                distributed_nparts = distributed ? nworkers : 1,
                distributed_partition_strategy = DISTRIBUTED_PARTITION_STRATEGY,
                threaded = USE_THREADED_PER_WORKER,
                state_filter = state_filter,
                state_input_filter = state_input_filter,
                print_level = PRINT_LEVEL,
            )
        end

        # 6) optimize!
        t_opt_wall = @elapsed MOI.optimize!(optimizer)

        t_construct = try
            MOI.get(optimizer, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
        catch
            NaN
        end

        global _last_benchmark_result = (
            nworkers = nworkers,
            t_master_packages = MASTER_STARTUP.t_master_packages,
            t_master_includes = MASTER_STARTUP.t_master_includes,
            t_addprocs = worker_setup.t_addprocs,
            t_worker_packages = worker_setup.t_worker_packages,
            t_worker_includes = worker_setup.t_worker_includes,
            t_warmup = t_warmup,
            t_problem_build = t_problem_build,
            t_optimizer_build = t_optimizer_build,
            t_opt_wall = t_opt_wall,
            t_abstract_reported = t_construct,
            t_end_to_end = 0.0,  # overwritten just below
            n_state = n_state,
            n_input = n_input,
        )
    end

    r = _last_benchmark_result
    return merge(r, (t_end_to_end = t_end_to_end,))
end

# ------------------------------------------------------------------------------
# Repeated benchmark
# ------------------------------------------------------------------------------

function benchmark_worker_count(nworkers::Int; nrepeat::Int = NREPEAT)
    results = NamedTuple[]

    println()
    println(
        "================================================================================",
    )
    println("Benchmarking with nworkers = $nworkers")
    println(
        "================================================================================",
    )

    for r in 1:nrepeat
        println("Run $r / $nrepeat")
        res = run_one_benchmark(nworkers)
        push!(results, res)

        @printf(
            "  addprocs = %.3f s | worker pkgs = %.3f s | worker includes = %.3f s | warmup = %.3f s\n",
            res.t_addprocs,
            res.t_worker_packages,
            res.t_worker_includes,
            res.t_warmup,
        )
        @printf(
            "  problem = %.3f s | optimizer = %.3f s | optimize! = %.3f s | reported abstraction = %.3f s | end-to-end = %.3f s\n",
            res.t_problem_build,
            res.t_optimizer_build,
            res.t_opt_wall,
            res.t_abstract_reported,
            res.t_end_to_end,
        )
    end

    return results
end

# ------------------------------------------------------------------------------
# Summary helpers
# ------------------------------------------------------------------------------

mean_or_nan(v) = isempty(v) ? NaN : mean(v)
std_or_nan(v) = length(v) <= 1 ? 0.0 : std(v)

function summarize_results(all_results::Vector{NamedTuple})
    worker_counts = sort(unique(r.nworkers for r in all_results))
    summary = NamedTuple[]

    for nw in worker_counts
        rows = filter(r -> r.nworkers == nw, all_results)

        push!(
            summary,
            (
                nworkers = nw,
                mean_addprocs = mean_or_nan([r.t_addprocs for r in rows]),
                std_addprocs = std_or_nan([r.t_addprocs for r in rows]),
                mean_worker_packages = mean_or_nan([r.t_worker_packages for r in rows]),
                std_worker_packages = std_or_nan([r.t_worker_packages for r in rows]),
                mean_worker_includes = mean_or_nan([r.t_worker_includes for r in rows]),
                std_worker_includes = std_or_nan([r.t_worker_includes for r in rows]),
                mean_warmup = mean_or_nan([r.t_warmup for r in rows]),
                std_warmup = std_or_nan([r.t_warmup for r in rows]),
                mean_problem_build = mean_or_nan([r.t_problem_build for r in rows]),
                std_problem_build = std_or_nan([r.t_problem_build for r in rows]),
                mean_optimizer_build = mean_or_nan([r.t_optimizer_build for r in rows]),
                std_optimizer_build = std_or_nan([r.t_optimizer_build for r in rows]),
                mean_opt_wall = mean_or_nan([r.t_opt_wall for r in rows]),
                std_opt_wall = std_or_nan([r.t_opt_wall for r in rows]),
                mean_abstract_reported = mean_or_nan([r.t_abstract_reported for r in rows]),
                std_abstract_reported = std_or_nan([r.t_abstract_reported for r in rows]),
                mean_end_to_end = mean_or_nan([r.t_end_to_end for r in rows]),
                std_end_to_end = std_or_nan([r.t_end_to_end for r in rows]),
            ),
        )
    end

    return summary
end

function print_summary(summary)
    println()
    println(
        "==============================================================================================================================",
    )
    println("Summary")
    println(
        "==============================================================================================================================",
    )
    @printf(
        "%7s | %18s | %18s | %18s | %18s | %18s\n",
        "workers",
        "addprocs mean±std",
        "worker pkgs mean±std",
        "worker inc mean±std",
        "warmup mean±std",
        "optimize! mean±std",
    )
    println("-"^130)

    for s in summary
        @printf(
            "%7d | %8.3f ± %-8.3f | %8.3f ± %-8.3f | %8.3f ± %-8.3f | %8.3f ± %-8.3f | %8.3f ± %-8.3f\n",
            s.nworkers,
            s.mean_addprocs,
            s.std_addprocs,
            s.mean_worker_packages,
            s.std_worker_packages,
            s.mean_worker_includes,
            s.std_worker_includes,
            s.mean_warmup,
            s.std_warmup,
            s.mean_opt_wall,
            s.std_opt_wall,
        )
    end

    println()
    @printf(
        "%7s | %18s | %18s | %18s | %18s\n",
        "workers",
        "problem mean±std",
        "optimizer mean±std",
        "reported abs mean±std",
        "end-to-end mean±std",
    )
    println("-"^105)

    for s in summary
        @printf(
            "%7d | %8.3f ± %-8.3f | %8.3f ± %-8.3f | %8.3f ± %-8.3f | %8.3f ± %-8.3f\n",
            s.nworkers,
            s.mean_problem_build,
            s.std_problem_build,
            s.mean_optimizer_build,
            s.std_optimizer_build,
            s.mean_abstract_reported,
            s.std_abstract_reported,
            s.mean_end_to_end,
            s.std_end_to_end,
        )
    end
end

function save_summary(summary, all_results; filename = RESULTS_FILE)
    open(filename, "w") do io
        println(io, "# Benchmark abstraction scaling")
        println(io, "# Generated: ", now())
        println(io, "# Master startup (measured once per Julia session)")
        println(io, "# t_master_packages = ", MASTER_STARTUP.t_master_packages)
        println(io, "# t_master_includes = ", MASTER_STARTUP.t_master_includes)
        println(io)

        println(io, "# Summary")
        for s in summary
            println(
                io,
                join(
                    [
                        string(s.nworkers),
                        string(s.mean_addprocs),
                        string(s.std_addprocs),
                        string(s.mean_worker_packages),
                        string(s.std_worker_packages),
                        string(s.mean_worker_includes),
                        string(s.std_worker_includes),
                        string(s.mean_warmup),
                        string(s.std_warmup),
                        string(s.mean_problem_build),
                        string(s.std_problem_build),
                        string(s.mean_optimizer_build),
                        string(s.std_optimizer_build),
                        string(s.mean_opt_wall),
                        string(s.std_opt_wall),
                        string(s.mean_abstract_reported),
                        string(s.std_abstract_reported),
                        string(s.mean_end_to_end),
                        string(s.std_end_to_end),
                    ],
                    ",",
                ),
            )
        end

        println(io)
        println(io, "# Raw results")
        for r in all_results
            println(
                io,
                join(
                    [
                        string(r.nworkers),
                        string(r.t_master_packages),
                        string(r.t_master_includes),
                        string(r.t_addprocs),
                        string(r.t_worker_packages),
                        string(r.t_worker_includes),
                        string(r.t_warmup),
                        string(r.t_problem_build),
                        string(r.t_optimizer_build),
                        string(r.t_opt_wall),
                        string(r.t_abstract_reported),
                        string(r.t_end_to_end),
                        string(r.n_state),
                        string(r.n_input),
                    ],
                    ",",
                ),
            )
        end
    end
end

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

function main()
    println("Master startup:")
    @printf("  master package load:  %.3f s\n", MASTER_STARTUP.t_master_packages)
    @printf("  master file includes: %.3f s\n", MASTER_STARTUP.t_master_includes)

    all_results = NamedTuple[]

    for nw in WORKER_COUNTS
        res = benchmark_worker_count(nw; nrepeat = NREPEAT)
        append!(all_results, res)
    end

    summary = summarize_results(all_results)
    print_summary(summary)

    if SAVE_RESULTS
        save_summary(summary, all_results)
        println("\nSaved results to: ", RESULTS_FILE)
    end

    reset_workers!()
    return summary, all_results
end

summary, all_results = main()
