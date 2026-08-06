module ThermostatHybridSystem

using StaticArrays
using LinearAlgebra
using Plots

import HybridSystems as HS
import MathematicalSystems as MS
import LazySets

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

Base.@kwdef struct Params{T}
    Ta::T = 18.0
    alpha::T = 0.1
    beta::T = 2.0
end

function off_dynamics(params::Params = Params())
    return (x, u) -> SVector(-params.alpha * (x[1] - params.Ta))
end

function on_dynamics(params::Params = Params())
    return (x, u) -> SVector(-params.alpha * (x[1] - params.Ta) + params.beta * u[1])
end

"""
    system(; params, _X_, _Uoff_, _Uon_)

Thermostat hybrid system whose modes are **plain physical systems**: the augmented
state is `(temperature, mode)`, and guards and the reset act on the temperature alone.
`ThermostatHybridTimeSystem` adds a clock subsystem (augmented state
`(temperature, time, mode)`).
"""
function system(;
    params::Params = Params(),
    _X_ = UT.box(SVector(17.0), SVector(25.0)),
    _Uoff_ = UT.box(SVector(0.0), SVector(0.0)),
    _Uon_ = UT.box(SVector(0.2), SVector(1.0)),
)
    off_system = MS.ConstrainedBlackBoxControlContinuousSystem(
        off_dynamics(params),
        1,
        1,
        _X_,
        _Uoff_,
    )
    on_system =
        MS.ConstrainedBlackBoxControlContinuousSystem(on_dynamics(params), 1, 1, _X_, _Uon_)

    # Plain modes — no `VectorContinuousSystem`, no clock subsystem.
    modes_systems = [off_system, on_system]

    reset_maps = [ST.GuardedResetMap(_X_), ST.GuardedResetMap(_X_)]  # identity on T

    automaton = HS.GraphAutomaton(2)
    HS.add_transition!(automaton, 1, 2, 1)
    HS.add_transition!(automaton, 2, 1, 2)
    switchings = [HS.AutonomousSwitching(), HS.AutonomousSwitching()]

    return HS.HybridSystem(automaton, modes_systems, reset_maps, switchings)
end

function jacobian_bounds(params::Params = Params())
    return [u -> SMatrix{1, 1}(-params.alpha), u -> SMatrix{1, 1}(-params.alpha)]
end

function make_cost_function(;
    switch_cost = 100.0,
    input_weight = 100.0,
    exp_weight = 100.0,
    exp_gain = 1000.0,
    u_soft = 0.4,
)
    # The cost depends only on the applied input `u`, not on the augmented state, so the
    # same function works whether or not the state carries a clock.
    return function (aug_state, u)
        if isa(u, String) && occursin("SWITCH", u)
            return switch_cost
        end

        uval = isa(u, Number) ? Float64(u) : isa(u, AbstractVector) ? Float64(u[1]) : 0.0

        energy_cost = input_weight * uval^2

        # Exponential penalty only above u_soft
        excess = max(uval - u_soft, 0.0)
        exp_cost = exp_weight * (exp(exp_gain * excess) - 1.0)

        return energy_cost + exp_cost
    end
end

"""
    optimal_control_problem(; params, initial_temperature, initial_mode, target)

Reach `target` (a temperature interval) in either mode over the [`system`](@ref) hybrid.
The augmented state is `(temperature, mode)` and the target is a mode-indexed `StateSpec`.
The time-windowed variant is `ThermostatHybridTimeSystem`.
"""
function optimal_control_problem(;
    params::Params = Params(),
    initial_temperature = 18.0,
    initial_mode = 1,
    target = UT.box(SVector(21.0), SVector(23.0)),
)
    hybrid_system = system(; params = params)

    initial_state = (SVector(initial_temperature), initial_mode)  # (temperature, mode)

    target_set = PR.HybridSpec(Dict(1 => PR.StateSpec(target), 2 => PR.StateSpec(target)))

    return PR.OptimalControlProblem(
        hybrid_system,
        initial_state,
        target_set,
        nothing,
        make_cost_function(),
    )
end

problem(; kwargs...) = optimal_control_problem(; kwargs...)

# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

include("thermostat_view.jl")

"""
    system_plot!(; problem, tlims = _THERMO_TLIMS)

Return a per-frame drawer `(fig, x, u, k) -> fig` for
[`Dionysos.animate_trajectory_dashboard`](@ref): the shared thermostat thermometer
(`thermostat_view.jl`), with the heater state read from the active mode `k`
(1 = OFF, 2 = ON) and the target band from the mode-indexed spec.
"""
function system_plot!(; problem, tlims = _THERMO_TLIMS)
    target = problem.target_set.per_mode[1].set
    tlo = LazySets.low(target, 1)
    thi = LazySets.high(target, 1)
    return function (fig, x, u, k)
        on = k == 2   # 1 = OFF, 2 = ON
        return _draw_thermometer!(
            fig,
            x[1],
            on,
            tlo,
            thi;
            tlims = tlims,
            label = "heater $(on ? "ON" : "OFF")",
        )
    end
end

end
