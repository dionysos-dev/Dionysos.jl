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

const TSTEP = 0.1
const DT_REF = 1e-4
const DT_TESTS = [5e-4, 1e-3, 5e-3]
const N_STEPS_TRAJ = 300

const X0 = SVector{6, Float64}(
    -0.10471975511965978,
    0.10471975511965978,
    0.10471975511965978,
    -0.45,
    0.0,
    0.0,
)

const U0 = SVector{3, Float64}(0.0, 0.0, 1.0)

function build_system(tstep, Δt_simu)
    problem = RobotProblem.problem(;
        robot_urdf = ROBOT_URDF_FILE,
        tstep = tstep,
        Δt_simu = Δt_simu,
        simulator = :custom, # :history,
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

function position_error(xs_ref, xs_test)
    n = min(length(xs_ref), length(xs_test))

    errs = Float64[]

    for k in 1:n
        # Compare only positions: first 3 states
        e = maximum(abs.(xs_ref[k][1:3] .- xs_test[k][1:3]))
        push!(errs, e)
    end

    return (max_error = maximum(errs), mean_error = mean(errs), final_error = errs[end])
end

function benchmark_accuracy()
    println("===============================================================")
    println("JuliaRobotics accuracy benchmark")
    println("===============================================================")
    println("Reference Δt_simu: ", DT_REF)
    println("TSTEP: ", TSTEP)
    println("N_STEPS_TRAJ: ", N_STEPS_TRAJ)
    println("X0: ", X0)
    println("U0: ", U0)
    println()

    println("Building reference system...")
    ref_system = build_system(TSTEP, DT_REF)

    println("Computing reference trajectory...")
    t_ref = @elapsed xs_ref = rollout(ref_system, X0, U0; nsteps = N_STEPS_TRAJ)
    @printf("Reference trajectory time: %.6f s\n", t_ref)
    println()

    for Δt_simu in DT_TESTS
        println("---------------------------------------------------------------")
        @printf("Testing Δt_simu = %.6g\n", Δt_simu)
        @printf("Internal steps per Dionysos step ≈ %.1f\n", TSTEP / Δt_simu)

        system = build_system(TSTEP, Δt_simu)

        t_test = @elapsed xs_test = rollout(system, X0, U0; nsteps = N_STEPS_TRAJ)

        err = position_error(xs_ref, xs_test)

        @printf("Trajectory time:       %.6f s\n", t_test)
        @printf("Speedup vs reference:  %.2fx\n", t_ref / t_test)
        @printf("Position max error:    %.8e rad\n", err.max_error)
        @printf("Position mean error:   %.8e rad\n", err.mean_error)
        @printf("Position final error:  %.8e rad\n", err.final_error)
    end

    println("===============================================================")
    println("Done")
    return println("===============================================================")
end

benchmark_accuracy()

# Δt_simu = 5e-4  → 5.7× faster, max error ≈ 0.013 rad ≈ 0.76° => good
# Δt_simu = 1e-3  → 10× faster, max error ≈ 0.299 rad ≈ 17.1°
# Δt_simu = 5e-3  → unstable / NaN
