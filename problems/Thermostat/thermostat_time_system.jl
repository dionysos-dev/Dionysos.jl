module ThermostatTimeSystem

# Continuous thermostat with a *time lift*: the state is augmented with a clock, so it is
# `(temperature, time)` and the dynamics are the switched thermostat plus `ṫ = 1`. This is
# the continuous analogue of `ThermostatHybridTimeSystem` (a clock added to the hybrid),
# completing the {continuous, hybrid} × {no lift, time lift} matrix. The input `u` selects
# the heater (`u = 1` ON, `u = 2` OFF), exactly as in `ThermostatSystem`.

import LazySets
using StaticArrays
using MathematicalSystems
using Plots
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

Base.@kwdef struct Params{T}
    Ta::T = 20.0
    alpha::T = 0.1
    beta::T = 2.0
end

_A(p::Params) = -p.alpha
_b_on(p::Params) = p.alpha * p.Ta + p.beta
_b_off(p::Params) = p.alpha * p.Ta

# Dynamics on the lifted state `(T, t)`: switched temperature dynamics plus `ṫ = 1`.
function dynamic(p::Params = Params())
    A, bon, boff = _A(p), _b_on(p), _b_off(p)
    return (x, u) -> SVector(u[1] == 1 ? A * x[1] + bon : A * x[1] + boff, 1.0)
end

# Jacobian of the dynamics w.r.t. `(T, t)`: `[-α 0; 0 0]` (constant, input-independent).
function jacobian_bound(p::Params = Params())
    J = SMatrix{2, 2}(-p.alpha, 0.0, 0.0, 0.0)
    return u -> J
end

function system(;
    params::Params = Params(),
    _X_ = UT.box(SVector(17.0, 0.0), SVector(25.0, 5.0)),   # (T, t)
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

"""
    optimal_control_problem(; params, _I_, _T_target, deadline, _X_, time)

Reach the temperature band `_T_target` while the clock is in `[min_time, max_time]`, over
the time-lifted `(temperature, time)` [`system`](@ref). The reach target is the box
`_T_target × [min_time, max_time]` in `(T, t)`, so the clock state enforces the timing:
`min_time = 0` is a plain deadline (reach before `max_time`), while a positive `min_time`
is a genuine window ("reach no earlier than `min_time`").

A positive `min_time` requires the abstraction's **time grid to be finer than the
integration step** (`h_t < Δt`); otherwise the GROWTH over-approximation lets the clock
"stall" and the window is (soundly but spuriously) unreachable — see the runner script.
"""
function optimal_control_problem(;
    params::Params = Params(),
    _I_ = UT.box(SVector(18.0, 0.0), SVector(18.5, 0.05)),
    _T_target = UT.box(SVector(21.0), SVector(23.0)),   # temperature band
    min_time = 3.0,
    max_time = 4.0,
    _X_ = UT.box(SVector(17.0, 0.0), SVector(25.0, 5.0)),
    time = DI.Problem.Infinity(),
)
    sys = system(; params = params, _X_ = _X_)

    # Target in `(T, t)`: reach the temperature band with the clock in `[min_time, max_time]`.
    target = UT.box(
        SVector(LazySets.low(_T_target, 1), min_time),
        SVector(LazySets.high(_T_target, 1), max_time),
    )

    state_cost = x -> 0.0
    transition_cost = (x, u) -> 1.0

    return PR.OptimalControlProblem(sys, _I_, target, state_cost, transition_cost, time)
end

problem(; kwargs...) = optimal_control_problem(; kwargs...)

# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

include("thermostat_view.jl")

"""
    system_plot!(; problem, tlims = _THERMO_TLIMS)

Return a per-frame drawer `(fig, x, u) -> fig` for
[`Dionysos.animate_trajectory_dashboard`](@ref): the shared thermostat thermometer
(`thermostat_view.jl`). The lifted state is `(temperature, time)`, so the temperature is
`x[1]`; the heater comes from the input (`u = 1` ON, `u = 2` OFF) and the target band
from the temperature range of `problem.target_set`.
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
