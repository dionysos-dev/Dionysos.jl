include(joinpath(@__DIR__, "..", "ControllerRuntime.jl"))
using .ControllerRuntime

cs = ControllerRuntime.Server.ControllerServer
(status, kind, version) =
    cs.load_controller!(; path = joinpath(@__DIR__, "..", "pid_cfg.jld2"))
println("Loaded controller: ", (status, kind, version))

for i in 1:5
    res = cs.step_controller(; x = [0.2, 0.0])
    println("Step $i: u = ", res.u)
end

(status,) = cs.reset_controller!()
println("Reset controller status: ", status)

for i in 1:5
    res = cs.step_controller(; x = [0.2, 0.0])
    println("Step $i: u = ", res.u)
end

(status,) = cs.unload_controller!()
println("Unloaded controller status: ", status)
