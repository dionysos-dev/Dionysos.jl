# Velocity-controlled pendulum: the joint is driven by a velocity servo, so the commanded angular
# velocity is the input and the model is the one-dimensional integrator `θ̇ = u` (N = M = 1, so
# the reduced model of the velocity-control reduction has dimension 2N − M = 1).
#
# Gravity has left the dynamics: it reappears only in the torque the servo must produce,
# `τ / (m l²) = θ̈ + (g/l) sin θ` (θ from the bottom, units of `SimplePendulum`, where the input is
# torque divided by m l²). A torque limit `U` therefore becomes two checks:
#
# - hold: between switches θ̈ = 0, so `|(g/l) sin θ| ≤ U` along the step. It depends on the state
#   and the current input only, and prunes transitions at abstraction time.
# - ramp: a switch u⁻ → u is commanded as a linear ramp of duration `t_r`, so
#   `|(u − u⁻)/t_r + (g/l) sin θ| ≤ U` during the ramp. It depends on the previous input and is
#   enforced by the pair synthesis.
#
# Both are checked over every angle the pendulum can visit, not only cell centers. With `t_r = τ`
# there is no hold phase, and holding a speed is the ramp check with `u = u⁻`.
module VelocityPendulum

import LazySets
using StaticArrays
using MathematicalSystems
using Dionysos
const PR = Dionysos.Problem
const MP = Dionysos.Mapping
const OPDS = Dionysos.Optim.DiscreteSystems

Base.@kwdef struct Params{T}
    l::T = 1.0
    g::T = 9.81
end

"""
    sin_range(lo, hi) -> (min, max)

Exact range of `sin` over `[lo, hi]`.
"""
function sin_range(lo::Real, hi::Real)
    hi - lo >= 2π && return (-1.0, 1.0)
    vals = [sin(lo), sin(hi)]
    k = ceil((lo - π / 2) / π)
    while π / 2 + k * π <= hi
        push!(vals, sin(π / 2 + k * π))
        k += 1
    end
    return extrema(vals)
end

# Range of the gravity term (g/l) sin θ over an interval of angles.
gravity_range(params::Params, lo, hi) = (params.g / params.l) .* sin_range(lo, hi)

"""
    discretization(; τ = 0.1, n = 125, levels = 4)

Commensurable grids on the circle: `n` cells of width `h = 2π/n`, and speed step `δu = h/τ`, so
that `θ̇ = u` maps cells exactly onto cells, including across the wrap at ±π. The cells tile
`[−π, π)` from −π; with `n` odd, θ = 0 is a cell center.
"""
function discretization(; τ = 0.1, n = 125, levels = 4)
    isodd(n) ||
        error("n must be odd, so that the cells tile [−π, π) with θ = 0 at a center")
    h = 2π / n
    δu = h / τ
    return (;
        τ,
        h,
        δu,
        u_max = levels * δu,
        # The periodic mapping wants the origin at the first cell center; with n odd this is the
        # same lattice as centers on multiples of h.
        state_grid = MP.GridFree(SVector(-π + h / 2), SVector(h)),
        input_grid = MP.GridFree(SVector(0.0), SVector(δu)),
    )
end

"""
    system(disc; θ_low = -π, θ_high = π)

The velocity-input pendulum `θ̇ = u` on the full circle; abstract it with a periodic mapping of
period 2π starting at −π, as the torque model.
"""
function system(disc; θ_low = -π, θ_high = π)
    _X_ = LazySets.Hyperrectangle(; low = SVector(θ_low), high = SVector(θ_high))
    _U_ = LazySets.Hyperrectangle(; low = SVector(-disc.u_max), high = SVector(disc.u_max))
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        (x, u) -> SVector(u[1]),
        1,
        1,
        _X_,
        _U_,
    )
end

"""
    hold_filter(params, disc, U)

`state_input_filter` for the hold phase: keep the transition `(θ, u)` only if the servo can hold
`u` against gravity at every angle swept during the step, from anywhere in the cell of `θ`.
"""
function hold_filter(params::Params, disc, U)
    return function (x, u)
        a, b = x[1], x[1] + disc.τ * u[1]
        lo, hi = gravity_range(params, min(a, b) - disc.h / 2, max(a, b) + disc.h / 2)
        return max(abs(lo), abs(hi)) <= U
    end
end

"""
    torque_limit(params, disc, U, t_r)

The ramp-phase torque limit as a state-dependent input constraint written on concrete cells and
inputs, for the grid abstraction to lift. The cell is that of the shifted angle `θ + (t_r/2) u⁻`
(the lattice is exact in this coordinate); during the ramp the true angle stays within
`(t_r/2) conv{−u⁻, 0, u}` of it, which is the exact range. Starts and ends at rest.
"""
function torque_limit(params::Params, disc, U, t_r)
    check = function (cell, u_prev, u)
        a = (u[1] - u_prev[1]) / t_r
        lo, hi = gravity_range(
            params,
            LazySets.low(cell)[1] + t_r / 2 * min(-u_prev[1], 0.0, u[1]),
            LazySets.high(cell)[1] + t_r / 2 * max(-u_prev[1], 0.0, u[1]),
        )
        return max(abs(a + lo), abs(a + hi)) <= U
    end
    return OPDS.StateInputVariation(
        check;
        target_input = SVector(0.0),
        initial_input = SVector(0.0),
    )
end

"""
    swing_up_problem(system; θ0 = 0.0, target_halfwidth = 15π/180, h = 0.05)

From rest at `θ0` (a set of width `h/2` around it, inside one cell of width `h`) to rest near
the top.
"""
function swing_up_problem(system; θ0 = 0.0, target_halfwidth = 15π / 180, h = 0.05)
    _I_ = LazySets.Hyperrectangle(; low = SVector(θ0 - h / 4), high = SVector(θ0 + h / 4))
    _T_ = LazySets.Hyperrectangle(;
        low = SVector(π - target_halfwidth),
        high = SVector(π + target_halfwidth),
    )
    return PR.OptimalControlProblem(system, _I_, _T_, nothing, (x, u) -> 1.0, PR.Infinity())
end

end # module
