include(joinpath(@__DIR__, "..", "ControllerRuntime.jl"))
include(joinpath(@__DIR__, "..", "transports", "lvserver_transport.jl"))

using .ControllerRuntime
using .LVServerTransport

# ControllerRuntime.Server.ControllerServer.load_controller!(;
#     path = joinpath(@__DIR__, "..", "pid_cfg.jld2"),
# )

LVServerTransport.run_lvserver(; initOK = true)
# LVServerTransport.run_json_server(endpoint="tcp://*:5556", verbose=true) # {"cmd":"load_controller","args":{"path":"C:/path/to/controller.file"}}, {"cmd":"step_controller","args":{"x":[0.1,0.2,0.3,0.4],"t":0.01}}
# LVServerTransport.run_both(initOK=true, json_endpoint="tcp://*:5556", verbose=true)

