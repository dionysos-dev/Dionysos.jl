using Printf
using Statistics
using Dionysos
using MathematicalSystems
using StaticArrays

const MS = MathematicalSystems

const ROBOT_PROBLEM_FILE = joinpath(@__DIR__, "..", "..", "6D_model", "robot_problem.jl")
const ROBOT_URDF_FILE =
    joinpath(@__DIR__, "..", "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf")

include(ROBOT_PROBLEM_FILE)
using .RobotProblem

# ==============================================================================
# Parameters
# ==============================================================================

const TSTEP = 0.1
const DT_REF = 1e-4
const DT_TESTS = [1e-4, 5e-4, 1e-3, 5e-3]
const BACKENDS = [:history, :custom]

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

# ==============================================================================
# Helpers
# ==============================================================================

function build_system(; tstep, Δt_simu, simulator)
    problem = RobotProblem.problem(;
        robot_urdf = ROBOT_URDF_FILE,
        tstep = tstep,
        Δt_simu = Δt_simu,
        simulator = simulator,
    )
    return problem.system
end

function rollout(system, x0, u; nsteps::Int)
    f = MS.mapping(system)

    xs = Vector{typeof(x0)}()
    sizehint!(xs, nsteps + 1)

    x = x0
    push!(xs, x)

    for _ in 1:nsteps
        x = f(x, u)
        push!(xs, x)
    end

    return xs
end

function first_bad_state(xs)
    for (k, x) in pairs(xs)
        if any(!isfinite, x)
            return k, x
        end
    end
    return nothing, nothing
end

function position_error(xs_ref, xs_test)
    n = min(length(xs_ref), length(xs_test))

    errs = Float64[]

    for k in 1:n
        e = maximum(abs.(xs_ref[k][1:3] .- xs_test[k][1:3]))
        push!(errs, e)
    end

    return (max_error = maximum(errs), mean_error = mean(errs), final_error = errs[end])
end

function benchmark_rollout(system, x0, u; nsteps::Int, nrepeats::Int)
    # Warm-up
    rollout(system, x0, u; nsteps = nsteps)

    times = Float64[]
    memories = Int[]

    for _ in 1:nrepeats
        t = @elapsed rollout(system, x0, u; nsteps = nsteps)
        mem = @allocated rollout(system, x0, u; nsteps = nsteps)

        push!(times, t)
        push!(memories, mem)
    end

    return (
        mean_time = mean(times),
        min_time = minimum(times),
        mean_memory = mean(memories),
        min_memory = minimum(memories),
    )
end

function print_row(backend, Δt_simu, stats, err, bad_step)
    @printf("%-8s  %.1e  ", string(backend), Δt_simu)

    if bad_step === nothing
        @printf(
            "%10.4f  %10.4f  %10.2f  %10.2f  %12.4e  %12.4e  %12.4e\n",
            stats.mean_time,
            stats.min_time,
            stats.mean_memory / 1024^2,
            stats.min_memory / 1024^2,
            err.max_error,
            err.mean_error,
            err.final_error,
        )
    else
        @printf(
            "%10.4f  %10.4f  %10.2f  %10.2f  %12s  %12s  %12s  bad@%s\n",
            stats.mean_time,
            stats.min_time,
            stats.mean_memory / 1024^2,
            stats.min_memory / 1024^2,
            "NaN",
            "NaN",
            "NaN",
            string(bad_step),
        )
    end
end

# ==============================================================================
# Main
# ==============================================================================

function run_comparison()
    println("===============================================================")
    println("History vs custom RigidBodyDynamics simulation benchmark")
    println("===============================================================")
    println("Reference: backend = :history, Δt_simu = ", DT_REF)
    println("TSTEP: ", TSTEP)
    println("N_STEPS_TRAJ: ", N_STEPS_TRAJ)
    println("N_REPEATS: ", N_REPEATS)
    println("X0: ", X0)
    println("U0: ", U0)
    println()

    println("Building reference system...")
    ref_system = build_system(; tstep = TSTEP, Δt_simu = DT_REF, simulator = :history)

    println("Computing reference trajectory...")
    ref_time = @elapsed xs_ref = rollout(ref_system, X0, U0; nsteps = N_STEPS_TRAJ)
    @printf("Reference rollout time: %.6f s\n", ref_time)
    println()

    println("Columns:")
    println(
        "backend   dt_simu      mean_time    min_time    mean_memMB   min_memMB      max_err      mean_err     final_err",
    )
    println(
        "---------------------------------------------------------------------------------------------------------------",
    )

    for backend in BACKENDS
        for Δt_simu in DT_TESTS
            system = build_system(; tstep = TSTEP, Δt_simu = Δt_simu, simulator = backend)

            stats = benchmark_rollout(
                system,
                X0,
                U0;
                nsteps = N_STEPS_TRAJ,
                nrepeats = N_REPEATS,
            )

            xs = rollout(system, X0, U0; nsteps = N_STEPS_TRAJ)

            bad_step, _ = first_bad_state(xs)

            if bad_step === nothing
                err = position_error(xs_ref, xs)
            else
                err = nothing
            end

            print_row(backend, Δt_simu, stats, err, bad_step)
        end
    end

    println("===============================================================")
    println("Done")
    return println("===============================================================")
end

run_comparison()

# Δt = 1e-4:
# speedup = 6.04 / 4.97 ≈ 1.21×
# memory reduction = 1082 / 1.5 ≈ 721×

# Δt = 5e-4:
# speedup = 1.09 / 0.99 ≈ 1.10×
# memory reduction = 234 / 1.5 ≈ 156×
