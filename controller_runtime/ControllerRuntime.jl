module ControllerRuntime

module Controllers
    include(joinpath(@__DIR__, "controllers", "artifacts.jl"))
    include(joinpath(@__DIR__, "controllers", "runtime_api.jl"))
    include(joinpath(@__DIR__, "controllers", "pid_artifact.jl"))
    include(joinpath(@__DIR__, "controllers", "pid_runtime.jl"))
end

module Server
    include(joinpath(@__DIR__, "server", "ControllerServer.jl"))
end

end # module
