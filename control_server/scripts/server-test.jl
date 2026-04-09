include(joinpath(@__DIR__, "..", "ControlServer.jl"))

using .ControlServer

server = ControlServer.ServerRuntime


f(x, u) = x[1] + u[1]
g(x) = x

controller = ControlServer.Controller4Server.Controller([(0.0)::Float64], f, g)

server.start_control_server(controller)