module PIDControllers

using LinearAlgebra
import Dionysos
const ST = Dionysos.System

export PIDMemory,
    PIDController,
    PIDControllerVector,
    ConstantSignal,
    ConstantTimeGetter,
    WrapAnglePositionVelocityError

struct PIDMemory{E, I}
    e_prev::E
    I::I
    initialized::Bool
end

# ------------------------------------------------------------------
# Small serializable helper callables
# ------------------------------------------------------------------

struct ConstantSignal{T}
    value::T
end
(c::ConstantSignal)(t) = c.value

struct ConstantTimeGetter end
(::ConstantTimeGetter)(x) = 0.0

struct WrapAnglePositionVelocityError end
function (e::WrapAnglePositionVelocityError)(x, r, t)
    eθ = mod(r[1] - x[1] + π, 2π) - π
    eω = r[2] - x[2]
    return typeof(r)(eθ, eω)
end

# ------------------------------------------------------------------
# PID controller
# ------------------------------------------------------------------

struct PIDController{K, EF, RF, DT, TG, U1, U2, E0} <: ST.AbstractContinuousController
    Kp::K
    Ki::K
    Kd::K
    error::EF
    ref::RF
    dt::DT
    time_getter::TG
    umin::U1
    umax::U2
    antiwindup::Bool
    e0::E0
end

_ref(ref, t) = ref
_ref(ref::ConstantSignal, t) = ref(t)

_dt(dt, x, t) = dt
_dt(dt::Function, x, t) = dt(x, t)

_apply_gain(K::AbstractMatrix, e) = K * e
_apply_gain(K::AbstractVector, e) = dot(K, e)
_apply_gain(K::Number, e) = K * e

_sat(u, umin, umax) = clamp.(u, umin, umax)

ST.controller_kind(::PIDController) = ST.DynamicKind()
ST.domain(pid::PIDController) = nothing

function ST.initial_state(pid::PIDController)
    return PIDMemory(pid.e0, zero(pid.e0), false)
end

function ST.output_control(pid::PIDController, mem::PIDMemory, x)
    t = pid.time_getter(x)
    r = _ref(pid.ref, t)
    e = pid.error(x, r, t)

    Δt = _dt(pid.dt, x, t)

    e_prev = mem.initialized ? mem.e_prev : e
    I_prev = mem.initialized ? mem.I : zero(e)

    I_new = I_prev + e * Δt
    de = (e - e_prev) / Δt

    u_unsat = _apply_gain(pid.Kp, e) + _apply_gain(pid.Ki, I_new) + _apply_gain(pid.Kd, de)

    if pid.umin !== nothing && pid.umax !== nothing
        return _sat(u_unsat, pid.umin, pid.umax)
    else
        return u_unsat
    end
end

function ST.update_state(pid::PIDController, mem::PIDMemory, x)
    t = pid.time_getter(x)
    r = _ref(pid.ref, t)
    e = pid.error(x, r, t)

    Δt = _dt(pid.dt, x, t)

    e_prev = mem.initialized ? mem.e_prev : e
    I_prev = mem.initialized ? mem.I : zero(e)

    I_new = I_prev + e * Δt
    de = (e - e_prev) / Δt

    u_unsat = _apply_gain(pid.Kp, e) + _apply_gain(pid.Ki, I_new) + _apply_gain(pid.Kd, de)

    if pid.umin !== nothing && pid.umax !== nothing
        u = _sat(u_unsat, pid.umin, pid.umax)
        if pid.antiwindup && any(u .!= u_unsat)
            I_new = I_prev
        end
    end

    return PIDMemory(e, I_new, true)
end

function PIDControllerVector(;
    Kp,
    Ki,
    Kd = zero(Kp),
    ref,
    error,
    dt = 1.0,
    time_getter = ConstantTimeGetter(),
    umin = nothing,
    umax = nothing,
    antiwindup::Bool = true,
    e0,
)
    return PIDController(
        Kp,
        Ki,
        Kd,
        error,
        ref,
        dt,
        time_getter,
        umin,
        umax,
        antiwindup,
        e0,
    )
end

end
