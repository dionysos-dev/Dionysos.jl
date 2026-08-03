module ThermostatSystem

import LazySets
using StaticArrays
using MathematicalSystems
using Plots
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

# ------------------------------------------------------------
# Model
# ------------------------------------------------------------
#
# State:
#   x = [T]
#   T = room temperature
#
# Input:
#   u = [1]  heater ON
#   u = [2]  heater OFF
#
# Continuous-time dynamics:
#
#   Ṫ = -α(T - Ta) + βq
#
# with:
#   Ta = ambient temperature
#   α  = heat-loss coefficient
#   β  = heater power
#   q  = 1 if heater ON, 0 if heater OFF
#
# Equivalently, this is a switched affine system:
#
#   heater ON:   Ṫ = -αT + αTa + β
#   heater OFF:  Ṫ = -αT + αTa
#
# This is analogous to the DCDC example: the input selects one
# continuous vector field.

Base.@kwdef struct Params{T}
    Ta::T = 20.0
    alpha::T = 0.1
    beta::T = 2.0
end

function A(p::Params = Params())
    return SMatrix{1, 1}(-p.alpha)
end

function b_on(p::Params = Params())
    return SVector(p.alpha * p.Ta + p.beta)
end

function b_off(p::Params = Params())
    return SVector(p.alpha * p.Ta)
end

function jacobian_bound(p::Params = Params())
    Ap = A(p)
    return u -> Ap
end

function DF_sys(p::Params = Params())
    Ap = A(p)
    return u -> Ap
end

function dynamic(p::Params = Params())
    Ap = A(p)
    bon = b_on(p)
    boff = b_off(p)

    return (x, u) -> u[1] == 1 ? Ap * x + bon : Ap * x + boff
end

function system(;
    params::Params = Params(),
    _X_ = UT.box(SVector(18.0), SVector(24.0)),
    _U_ = UT.box(SVector(1), SVector(2)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        LazySets.dim(_X_),
        LazySets.dim(_U_),
        _X_,
        _U_,
    )
end

# ------------------------------------------------------------
# Problems
# ------------------------------------------------------------

"""
    safety_problem(; params, _I_, _S_, time)

Keep the room temperature inside the safe interval `_S_`.

Default:
- initial temperature: 20°C to 21°C
- safe temperature: 18°C to 24°C
"""
function safety_problem(;
    params::Params = Params(),
    _I_ = UT.box(SVector(20.0), SVector(21.0)),
    _S_ = UT.box(SVector(18.0), SVector(24.0)),
    time = PR.Infinity(),
)
    sys = system(; params = params, _X_ = _S_)

    return PR.SafetyProblem(sys, _I_, _S_, time)
end

"""
    optimal_control_problem(; params, _I_, _T_, _X_, time)

Finite-horizon reachability problem.

Reach the target temperature interval `_T_` from `_I_` in at most `time`.
The costs are set to zero, so this is pure reachability.
"""
function optimal_control_problem(;
    params::Params = Params(),
    _I_ = UT.box(SVector(18.0), SVector(18.5)),
    _T_ = UT.box(SVector(21.8), SVector(22.2)),
    _X_ = UT.box(SVector(18.0), SVector(24.0)),
    time = 30.0,
)
    sys = system(; params = params, _X_ = _X_)

    state_cost = x -> 0.0
    transition_cost = (x, u) -> 1.0

    return PR.OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost, time)
end

"""
    reach_and_stay_problem(; params, _I_, _T_, _S_, time)

Reach the target temperature interval `_T_` and stay there forever,
while remaining inside the safe set `_S_`.
"""
function reach_and_stay_problem(;
    params::Params = Params(),
    _I_ = UT.box(SVector(18.0), SVector(18.5)),
    _T_ = UT.box(SVector(21.8), SVector(22.2)),
    _S_ = UT.box(SVector(18.0), SVector(24.0)),
    time = PR.Infinity(),
)
    sys = system(; params = params, _X_ = _S_)

    return PR.ReachAndStayProblem(sys, _I_, _T_, _S_, time)
end

# Keep the same naming convention as your DCDC module.
problem(; kwargs...) = safety_problem(; kwargs...)

# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

include("thermostat_view.jl")

"""
    system_plot!(; problem, tlims = _THERMO_TLIMS)

Return a per-frame drawer `(fig, x, u) -> fig` for
[`Dionysos.animate_trajectory_dashboard`](@ref): the shared thermostat thermometer
(`thermostat_view.jl`). The continuous input selects the heater (`u = 1` ON, `u = 2`
OFF); the target band is read from `problem.target_set` (a reach / reach-and-stay
problem).
"""
function system_plot!(; problem, tlims = _THERMO_TLIMS)
    tlo = LazySets.low(problem.target_set, 1)
    thi = LazySets.high(problem.target_set, 1)
    return function (fig, x, u)
        on = Int(round(u[1])) == 1   # continuous input: 1 = ON, 2 = OFF
        return _draw_thermometer!(
            fig,
            Float64(x[1]),
            on,
            tlo,
            thi;
            tlims = tlims,
            label = "heater $(on ? "ON" : "OFF")",
        )
    end
end

end
