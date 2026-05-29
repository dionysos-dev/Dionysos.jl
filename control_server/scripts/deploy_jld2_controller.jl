include(joinpath(@__DIR__, "..", "ControlServer.jl"))

using .ControlServer
using JLD2

const CSD = ControlServer.ControlServerDeployment
const SR = ControlServer.ServerRuntime

# ------------------------------------------------------------
# Path to saved controller
# ------------------------------------------------------------
# controller_file = joinpath(@__DIR__, "pid_controller.jld2")

include(
    joinpath(
        @__DIR__,
        "..",
        "..",
        "problems",
        "BipedRobot",
        "6D_model",
        "two_step_walking_controller.jl",
    ),
)
controller_file = joinpath(@__DIR__, "robot_two_step_controller_CM.jld2")

# ------------------------------------------------------------
# Load Dionysos controller and convert for server runtime
# ------------------------------------------------------------
server_controller = CSD.load_server_controller(controller_file)

# ------------------------------------------------------------
# Start server
# ------------------------------------------------------------
port = 5000

result = SR.start_control_server(
    server_controller;
    port = port,
    log_data = true,
    received_data_size = 8,
    state_to_vector = nothing,
)

# ------------------------------------------------------------
# Optional post-processing
# ------------------------------------------------------------
if result !== nothing
    t, measurements, controls, states = result
    println("Server session finished.")
    println("Logged $(length(t)) packets.")
    println("Average dt: $((t[length(t)]-t[1])/length(t))")
end

using Plots

p = plot(; layout = (2, 2), size = (1200, 900))

for i in 1:4
    plot!(p[i], t, measurements[i, :]; label = "Angle $i")

    #= plot!(
        p[i],
        t,
        measurements[i+4, :],
        label = "Velocity $i"
    ) =#

    plot!(p[i], t, controls[i, :]; label = "Control input $i")

    xlabel!(p[i], "Time")
end

display(p)
