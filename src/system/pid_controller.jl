module PIDControllers
using LinearAlgebra
import MathematicalSystems
const MS = MathematicalSystems

mutable struct PIDState
    e_prev::Any
    I::Any
    initialized::Bool
end

# K can be scalar, vector, or matrix
mutable struct PIDController{K}
    Kp::K
    Ki::K
    Kd::K
    error::Any        # (x, r, t) -> e   (vector)
    ref::Any          # constant ref or ref(t)
    dt::Any           # constant dt or dt(x,t)
    umin::Any
    umax::Any
    antiwindup::Bool
    st::PIDState
end

_ref(ref, t) = ref
_ref(ref::Function, t) = ref(t)

_dt(dt, x, t) = dt
_dt(dt::Function, x, t) = dt(x, t)

# multiply gain (scalar/vector/matrix) by vector error
_apply_gain(K::AbstractMatrix, e) = K * e     # m×n × n → m
_apply_gain(K::AbstractVector, e) = dot(K, e) # 1×n × n → scalar
_apply_gain(K::Number, e) = K * e             # scalar error case

_sat(u, umin, umax) = clamp.(u, umin, umax)

function pid_map(
    pid::PIDController;
    nx::Int,
    nu::Int,
    time_getter = x -> 0.0,
    silent = true,
    X = nothing,
)
    h = function (x)
        t = time_getter(x)
        r = _ref(pid.ref, t)
        e = pid.error(x, r, t)

        dt = _dt(pid.dt, x, t)

        if !pid.st.initialized
            pid.st.e_prev = e
            pid.st.I = zero(e)
            pid.st.initialized = true
        end

        pid.st.I = pid.st.I + e * dt
        de = (e - pid.st.e_prev) / dt

        u_unsat =
            _apply_gain(pid.Kp, e) + _apply_gain(pid.Ki, pid.st.I) + _apply_gain(pid.Kd, de)

        u = u_unsat
        if pid.umin !== nothing && pid.umax !== nothing
            u = _sat(u_unsat, pid.umin, pid.umax)
            if pid.antiwindup && any(u .!= u_unsat)
                pid.st.I = pid.st.I - e * dt
            end
        end

        if !silent
            println("e: ", e, " u: ", u)
            # println("i: ",  pid.st.I)
        end

        pid.st.e_prev = e
        return u
    end

    return MS.ConstrainedBlackBoxMap(nx, nu, h, X)
end

function PIDControllerVector(;
    Kp,
    Ki,
    Kd = zero(Kp),
    ref,
    error,
    dt = 1.0,
    umin = nothing,
    umax = nothing,
    antiwindup::Bool = true,
    e0,
)
    st = PIDState(e0, zero(e0), false)
    return PIDController(Kp, Ki, Kd, error, ref, dt, umin, umax, antiwindup, st)
end

end # module
