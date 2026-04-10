include(joinpath(@__DIR__, "..", "ControlServer.jl"))

using .ControlServer
using Plots

server = ControlServer.ServerRuntime

function f(x_c, x_p; Ts = 0.01::Float64)
    #     e  = -x_p[1]
    # int{e} = x_c

    return [x_c[1] - Ts*x_p[1]]
end

function g(x_c, x_p)
    #      e  = -x_p[1]
    # \dot{e} = -x_p[2]
    # \int{e} = x_c[1]

    kp = 0.6649::Float64
    ki = 0.2817::Float64
    kd = 0.3923::Float64

    return [-kp*x_p[1] + ki*x_c[1] - kd*x_p[2]]
end

controller = ControlServer.Controller4Server.Controller([(0.0)::Float64], f, g)

plotting = true

if plotting
    (t, xp, xc) = server.start_control_server(controller; log_data=true, received_data_size=2)
    println("Data stored")
    plot(t, [xp[1,:] xc[1,:]], layout=(2,1), label=["y", "u"])
else
    server.start_control_server(controller; log_data=false)
    println("Data not stored")
end