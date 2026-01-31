include(joinpath(@__DIR__, "..", "ControllerRuntime.jl"))
include(joinpath(@__DIR__, "..", "transports", "lvserver_transport.jl"))

using .ControllerRuntime
using .LVServerTransport

ControllerRuntime.Server.ControllerServer.load_controller!(;
    path = joinpath(@__DIR__, "..", "pid_cfg.jld2"),
)

LVServerTransport.start_lvserver(; initOK = true)
