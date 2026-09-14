using Dionysos
using StaticArrays
using JLD2
using LazySets

ST = Dionysos.System

include("../robot_vcontrol.jl")
import .RobotVelocity as RV

include("../../../../control_server/ControlServer.jl")
import .ControlServer

server = ControlServer.ServerRuntime

if !@isdefined(concrete_controller)
    # const concrete_controller = JLD2.load("./4D_model_vcontrol/controller.jld2", "controller")
    const concrete_controller = JLD2.load("./4D_model_vcontrol/experiment/pid_controller.jld2")
    println("Controller loaded from JLD2 file")
else
    println("Controller already in memory, skipping reload")
end

f(x_c, x_p) = [0.0]

function make_g(concrete_controller)

    function g(x_c, x_p)
        x_p_rad = x_p*pi/180
        u_rad = ST.output_control(concrete_controller, nothing, x_p_rad)
        return u_rad === nothing ? nothing : u_rad/pi*180
    end

    return g
end

controller = ControlServer.Controller4Server.Controller([0.0], f, make_g(concrete_controller))

result = server.start_control_server(
        controller;
        log_data = true,
        received_data_size = 4,
        state_to_vector = x -> x,
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
