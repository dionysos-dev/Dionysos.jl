include(joinpath(@__DIR__, "..", "ControllerRuntime.jl"))
using .ControllerRuntime
using StaticArrays

art = ControllerRuntime.Controllers.PIDArtifact.make_pid_artifact(;
    tstep = 0.1,
    Kp = SVector(15.0, 5.0),
    Ki = SVector(3.0, 0.0),
    Kd = SVector(0.0, 0.0),
    ref = SVector(3*pi/4, 0.0),
    umin = nothing, # 10.0
    umax = nothing, # 10.0
    antiwindup = true,
)

ControllerRuntime.Controllers.Artifacts.save_artifact(
    joinpath(@__DIR__, "..", "pid_cfg.jld2"),
    art,
)

println("Saved pid_cfg.jld2")
