using Printf
using Dionysos
using StaticArrays
using JLD2

const ST = Dionysos.System

const NSTEP = 300
const TSTEP = 0.1

const USE_TEST_TRAJECTORY = false
const VISUALIZE_TRAJECTORY = true

const CONTROLLER_FILE = joinpath(@__DIR__, "..", "..", "out", "6D", "robot_controller.jld2")
const ROBOT_URDF_FILE = joinpath(@__DIR__, "..", "..", "deps", "ZMP_2DBipedRobot_nodamping.urdf")
const RSVIZ_FILE = joinpath(@__DIR__, "..", "..", "src", "RSVisualization.jl")
include(RSVIZ_FILE)
using .RSVisualization

function load_controller_data(filename::AbstractString)
    return JLD2.jldopen(filename, "r") do file
        (
            controller = file["controller"],
            concrete_system = file["concrete_system"],
            control_problem = file["control_problem"],
            tstep = file["tstep"],
        )
    end
end

reached_target(problem) = x -> (x ∈ problem.target_set)

function make_test_trajectory(; N = 300, dt = 0.1)
    seq = Vector{SVector{6, Float64}}(undef, N)

    for k in 1:N
        t = (k - 1) * dt

        q3 = 5π / 180 * sin(0.5 * t)
        q4 = 4π / 180 * sin(0.5 * t + π / 4)
        q5 = 6π / 180 * sin(0.5 * t + π / 2)

        v3 = 5π / 180 * 0.5 * cos(0.5 * t)
        v4 = 4π / 180 * 0.5 * cos(0.5 * t + π / 4)
        v5 = 6π / 180 * 0.5 * cos(0.5 * t + π / 2)

        seq[k] = @SVector [q3, q4, q5, v3, v4, v5]
    end

    return ST.Trajectory(seq)
end

function visualize_trajectory!(x_traj; robot_urdf::AbstractString, tstep::Float64)
    rs, vis = RSVisualization.get_visualization_tool(; robot_urdf = robot_urdf)
    RSVisualization.animate_trajectory!(vis, x_traj.seq; dt = tstep)

    println("MeshCat visualizer opened. Keep this Julia session alive to view the animation.")
    return rs, vis
end

if USE_TEST_TRAJECTORY
    x_traj = make_test_trajectory(; N = NSTEP, dt = TSTEP)
    println("Synthetic trajectory length: ", length(x_traj.seq))
else
    @info "Loading controller" CONTROLLER_FILE
    data = load_controller_data(CONTROLLER_FILE)

    controller = data.controller
    concrete_system = data.concrete_system
    control_problem = data.control_problem
    tstep = data.tstep

    x0 = SVector{6, Float64}(zeros(6))

    @info "Simulating closed loop" NSTEP

    t_sim = @elapsed begin
        global x_traj, u_traj = ST.get_closed_loop_trajectory(
            concrete_system,
            controller,
            x0,
            NSTEP;
            stopping = reached_target(control_problem), # stopping = x -> false,
            verbose = true,
        )
    end

    @printf("Simulation time: %.3f s\n", t_sim)
end

println("Trajectory length: ", length(x_traj.seq))
println("Final state:")
println(last(x_traj.seq))

if VISUALIZE_TRAJECTORY
    global rs, vis = visualize_trajectory!(
        x_traj;
        robot_urdf = ROBOT_URDF_FILE,
        tstep = tstep,
    )
end

println("Done.")