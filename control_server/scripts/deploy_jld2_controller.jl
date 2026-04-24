include(joinpath(@__DIR__, "..", "ControlServer.jl"))

using .ControlServer
using JLD2

const CSD = ControlServer.ControlServerDeployment
const SR = ControlServer.ServerRuntime

# ------------------------------------------------------------
# Path to saved controller
# ------------------------------------------------------------
controller_file = joinpath(@__DIR__, "pid_controller.jld2")

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
    received_data_size = 2,
    state_to_vector = nothing,
)

# ------------------------------------------------------------
# Optional post-processing
# ------------------------------------------------------------
if result !== nothing
    t, measurements, controls, states = result
    println("Server session finished.")
    println("Logged $(length(t)) packets.")
end

using Plots

plot(t, measurements[1, :]; label = "y")
plot!(t, controls[1, :]; label = "u")
