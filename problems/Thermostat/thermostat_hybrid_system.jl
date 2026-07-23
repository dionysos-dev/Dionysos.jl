module ThermostatHybridSystem

using StaticArrays
using LinearAlgebra

import HybridSystems as HS
import MathematicalSystems as MS

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

struct ThermostatIdentityResetMap <: MS.AbstractMap
    domain::UT.Box
end

MS.apply(reset::ThermostatIdentityResetMap, state::AbstractVector) = state
MS.stateset(reset::ThermostatIdentityResetMap) = reset.domain

function off_dynamics(params::Params = Params())
    return (x, u) -> SVector(-params.alpha * (x[1] - params.Ta))
end

function on_dynamics(params::Params = Params())
    return (x, u) -> SVector(-params.alpha * (x[1] - params.Ta) + params.beta * u[1])
end

function system(;
    params::Params = Params(),
    _X_ = UT.box(SVector(17.0), SVector(25.0)),
    _Uoff_ = UT.box(SVector(0.0), SVector(0.0)),
    _Uon_ = UT.box(SVector(0.2), SVector(1.0)),
    _T_ = UT.box(SVector(0.0), SVector(100.0)),
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

    off_time_system = MS.ConstrainedLinearContinuousSystem(SMatrix{1, 1}(1.0), _T_)

    on_time_system = MS.ConstrainedLinearContinuousSystem(SMatrix{1, 1}(1.0), _T_)

    modes_systems = [
        ST.VectorContinuousSystem([off_system, off_time_system]),
        ST.VectorContinuousSystem([on_system, on_time_system]),
    ]

    switch_domain = UT.box(SVector(17.0, 0.0), SVector(25.0, 100.0))

    reset_maps = [
        ThermostatIdentityResetMap(switch_domain),
        ThermostatIdentityResetMap(switch_domain),
    ]

    automaton = HS.GraphAutomaton(2)

    HS.add_transition!(automaton, 1, 2, 1)
    HS.add_transition!(automaton, 2, 1, 2)

    switchings = [HS.AutonomousSwitching(), HS.AutonomousSwitching()]

    return HS.HybridSystem(automaton, modes_systems, reset_maps, switchings)
end

"""
    system_time_free(; params, _X_, _Uoff_, _Uon_)

Thermostat hybrid system whose modes are **plain physical systems** (no time
subsystem): the abstraction therefore carries no time axis and its augmented state
is `(temperature, mode)`. Guards and the reset act on the temperature alone.
"""
function system_time_free(;
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

    reset_maps = [ThermostatIdentityResetMap(_X_), ThermostatIdentityResetMap(_X_)]  # identity on T

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
    # `aug_state` is `(x, t, mode)` (timed) or `(x, mode)` (time-free); the cost only
    # depends on the applied input `u`, so it is arity-agnostic.
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

function optimal_control_problem(;
    params::Params = Params(),
    initial_temperature = 21.5,
    initial_mode = 1,
    target = UT.box(SVector(21.0), SVector(23.0)),
    target_time = UT.box(SVector(15.0), SVector(50.0)),
    time_domain = UT.box(SVector(0.0), SVector(50.0)),
)
    hybrid_system = system(; params = params, _T_ = time_domain)

    initial_state = (SVector(initial_temperature), 0.0, initial_mode)

    target_set = PR.hybrid_reach_spec([target, target], [target_time, target_time], [1, 2])

    transition_cost = make_cost_function(;
    # switch_cost = 1.0,
    # input_weight = 10000.0,
    # early_energy_weight = 0.0,
    # desired_time = 8.0,
    # u_soft = 0.2,
    # excess_weight = 1000,
    )

    return PR.OptimalControlProblem(
        hybrid_system,
        initial_state,
        target_set,
        nothing,
        transition_cost,
    )
end

problem(; kwargs...) = optimal_control_problem(; kwargs...)

"""
    optimal_control_problem_time_free(; params, initial_temperature, initial_mode, target)

Reach `target` (a temperature interval) in either mode, with **no time constraint**,
over the time-free [`system_time_free`](@ref) hybrid. The augmented state is
`(temperature, mode)` and the target is a mode-indexed [`StateSpec`](@ref) — no time
window is involved.
"""
function optimal_control_problem_time_free(;
    params::Params = Params(),
    initial_temperature = 18.0,
    initial_mode = 1,
    target = UT.box(SVector(21.0), SVector(23.0)),
)
    hybrid_system = system_time_free(; params = params)

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

problem_time_free(; kwargs...) = optimal_control_problem_time_free(; kwargs...)

end
