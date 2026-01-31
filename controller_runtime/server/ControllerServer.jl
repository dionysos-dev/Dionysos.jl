module ControllerServer

# Absolute references (robust during include)
import ...Controllers.Artifacts as CA
import ...Controllers.RuntimeAPI as RT
import ...Controllers.PIDRuntime as PID

export load_controller!, unload_controller!, reset_controller!, step_controller!

const _runtime = Ref{Union{Nothing, RT.AbstractRuntimeController}}(nothing)
const _loaded_path = Ref{Union{Nothing, String}}(nothing)

function load_controller!(; path::AbstractString)
    art = CA.load_artifact(path)

    rt = if art.kind == :pid
        PID.build_runtime(art)
    else
        error("Unsupported controller kind: $(art.kind)")
    end

    _runtime[] = rt
    _loaded_path[] = String(path)
    return (; ok = true, kind = String(art.kind), version = art.version)
end

function reset_controller!(; kwargs...)
    rt = _runtime[]
    rt === nothing && return (; ok = false, msg = "not loaded")
    return RT.reset!(rt; kwargs...)
end

function step_controller(; x, t = 0.0, kwargs...)
    rt = _runtime[]
    rt === nothing && return (; ok = false, msg = "not loaded")
    return RT.step!(rt; x = x, t = t, kwargs...)
end

function unload_controller!()
    _runtime[] = nothing
    _loaded_path[] = nothing
    return (; ok = true)
end

end # module
