module PIDRuntime

using StaticArrays
import Dionysos
const ST = Dionysos.System
const PID = ST.PIDControllers

import ..Artifacts: Artifact
import ..RuntimeAPI: AbstractRuntimeController, step!, reset!

export PIDRuntimeController, build_runtime

mutable struct PIDRuntimeController <: AbstractRuntimeController
    pid_obj::Any
    pid_map::Any
end

wrap_angle(e) = mod(e + π, 2π) - π
_to_svec2(x) = x isa SVector{2} ? x : SVector{2, Float64}(x)

function build_error()
    return (x, r, _) -> begin
        eθ = wrap_angle(r[1] - x[1])
        eω = r[2] - x[2]
        SVector(eθ, eω)
    end
end

function build_runtime(art::Artifact)::PIDRuntimeController
    p = art.payload

    antiwindup = get(p, "antiwindup", true)
    umin = get(p, "umin", nothing)
    umax = get(p, "umax", nothing)

    pid = PID.PIDControllerVector(;
        Kp = p["Kp"],
        Ki = p["Ki"],
        Kd = p["Kd"],
        ref = p["ref"],
        error = build_error(),
        dt = p["tstep"],
        umin = (umin === nothing ? nothing : SVector(Float64(umin))),
        umax = (umax === nothing ? nothing : SVector(Float64(umax))),
        antiwindup = antiwindup,
        e0 = p["e0"],
    )

    pid.st.I = p["I0"]
    pid.st.initialized = false

    pid_map = PID.pid_map(pid; nx = 2, nu = 1, silent = true)
    return PIDRuntimeController(pid, pid_map)
end

function step!(ctrl::PIDRuntimeController; x, t = 0.0, ref = nothing)
    xs = _to_svec2(x)
    if ref !== nothing
        ctrl.pid_obj.ref = _to_svec2(ref)
    end
    u = ctrl.pid_map.h(xs)
    return (; ok = true, u = u)
end

function reset!(ctrl::PIDRuntimeController; e0 = SVector(0.0, 0.0))
    ctrl.pid_obj.st.e_prev = e0
    ctrl.pid_obj.st.I = zero(e0)
    ctrl.pid_obj.st.initialized = false
    return (; ok = true)
end

end # module
