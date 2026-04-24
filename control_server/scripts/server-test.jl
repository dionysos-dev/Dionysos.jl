include(joinpath(@__DIR__, "..", "ControlServer.jl"))

using .ControlServer
using Plots

server = ControlServer.ServerRuntime

function f(x_c, x_p; Ts = 0.01)
    return [x_c[1] - Ts * x_p[1]]
end

function g(x_c, x_p)
    kp = 0.6649
    ki = 0.2817
    kd = 0.3923

    return [-kp * x_p[1] + ki * x_c[1] - kd * x_p[2]]
end

controller = ControlServer.Controller4Server.Controller([0.0], f, g)

plotting = true

if plotting
    t, xp, u, xc = server.start_control_server(
        controller;
        log_data = true,
        received_data_size = 2,
        state_to_vector = x -> x,
    )

    println("Data stored")

    plt = plot(; layout = (3, 1))
    plot!(plt[1], t, xp[1, :]; label = "y")
    plot!(plt[2], t, u[1, :]; label = "u")
    plot!(plt[3], t, xc[1, :]; label = "x_c")
    display(plt)
else
    server.start_control_server(controller; log_data = false)
    println("Data not stored")
end
