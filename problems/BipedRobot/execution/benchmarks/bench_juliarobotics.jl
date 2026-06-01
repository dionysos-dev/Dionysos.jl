using Printf
using Statistics
using Dionysos
using MathematicalSystems
using StaticArrays
using BenchmarkTools

const MS = MathematicalSystems
const ST = Dionysos.System

const ROBOT_PROBLEM_FILE = joinpath(@__DIR__, "..", "..", "6D_model", "robot_problem.jl")
const ROBOT_URDF_FILE =
    joinpath(@__DIR__, "..", "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf")

include(ROBOT_PROBLEM_FILE)
using .RobotProblem

const TSTEPS = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1]
const DT_SIMUS = [1e-4, 5e-4, 1e-3, 5e-3]
const N_STEPS_TRAJ = 300
const N_REPEATS = 5

const X0 = SVector{6, Float64}(
    -0.10471975511965978,
    0.10471975511965978,
    0.10471975511965978,
    -0.45,
    0.0,
    0.0,
)

const U0 = SVector{3, Float64}(0.0, 0.0, 1.0)

function build_problem(tstep, Δt_simu)
    return RobotProblem.problem(;
        robot_urdf = ROBOT_URDF_FILE,
        tstep = tstep,
        Δt_simu = Δt_simu,
        simulator = :custom, # :history,
    )
end

function benchmark_one_step(system, x, u)
    f = MS.mapping(system)

    f(x, u)

    b = @benchmark $f($x, $u) samples = 100 evals = 1

    return (
        median_time = median(b).time / 1e6,
        min_time = minimum(b).time / 1e6,
        allocations = median(b).allocs,
        memory = median(b).memory,
    )
end

function benchmark_trajectory(system, x0, u)
    f = MS.mapping(system)

    function rollout()
        x = x0
        for _ in 1:N_STEPS_TRAJ
            x = f(x, u)
        end
        return x
    end

    rollout()

    times = Float64[]
    allocs = Int[]

    for _ in 1:N_REPEATS
        t = @elapsed rollout()
        a = @allocated rollout()
        push!(times, t)
        push!(allocs, a)
    end

    return (
        mean_time = mean(times),
        min_time = minimum(times),
        allocations_bytes = minimum(allocs),
    )
end

function run_benchmark()
    println("===============================================================")
    println("JuliaRobotics / RobotProblem benchmark")
    println("===============================================================")
    println("N_STEPS_TRAJ: ", N_STEPS_TRAJ)
    println("X0: ", X0)
    println("U0: ", U0)
    println()

    for Δt_simu in DT_SIMUS
        println("===============================================================")
        @printf("Δt_simu = %.6g\n", Δt_simu)
        println("===============================================================")

        for tstep in TSTEPS
            ratio = tstep / Δt_simu

            if ratio < 1 || !isapprox(ratio, round(ratio); rtol = 1e-10, atol = 1e-12)
                println("---------------------------------------------------------------")
                @printf("tstep = %.6g\n", tstep)
                @printf("internal steps ≈ %.3f\n", ratio)
                println("SKIP: tstep must be an integer multiple of Δt_simu")
                continue
            end
            println("---------------------------------------------------------------")
            @printf("tstep = %.6g\n", tstep)
            @printf("internal steps ≈ %.1f\n", tstep / Δt_simu)

            problem = build_problem(tstep, Δt_simu)
            system = problem.system

            RobotProblem.warmup_robot_problem!(;
                robot_urdf = ROBOT_URDF_FILE,
                tstep = tstep,
                Δt_simu = Δt_simu,
            )

            one = benchmark_one_step(system, X0, U0)
            traj = benchmark_trajectory(system, X0, U0)

            @printf("one-step median:     %.6f ms\n", one.median_time)
            @printf("one-step min:        %.6f ms\n", one.min_time)
            @printf("one-step allocs:     %d\n", one.allocations)
            @printf("one-step memory:     %.3f KiB\n", one.memory / 1024)

            @printf("trajectory mean:     %.6f s\n", traj.mean_time)
            @printf("trajectory min:      %.6f s\n", traj.min_time)
            @printf("trajectory memory:   %.3f MiB\n", traj.allocations_bytes / 1024^2)
        end
    end

    println("===============================================================")
    println("Done")
    return println("===============================================================")
end

run_benchmark()
