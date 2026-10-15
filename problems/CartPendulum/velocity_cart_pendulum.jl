# Velocity-controlled cart-pole: a velocity servo drives the cart, so the cart velocity v is the
# input. With the pole angle θ, measured from the upright, and its conjugate momentum
# p = m l² θ̇ + m l cos θ v, the reduced model has state (x, θ, p), of dimension 2N − M = 3. The
# pole part does not depend on the cart position x, so without track limits (θ, p) is the whole
# state.
#
# Point-mass pole, no friction (Assumption 1 of the paper):
#   θ̇ = (p − m l cos θ v) / (m l²)
#   ṗ = −(sin θ / l) v p + m sin θ cos θ v² + m g l sin θ
# and the cart force that realizes a velocity ramp of slope a = v̇ is
#   F = (M + m sin² θ) a + m sin θ (g cos θ − l θ̇²).
#
# The servo ramps each command over the whole step, so the input during a step is a ramp from
# the previous command, not a constant. The abstraction is built for constant inputs, and the
# ramp enters as an additive disturbance |f(z, w) − f(z, u)| ≤ b |w − u| ≤ b Δ_max, with
# b = sup |∂f/∂w| and Δ_max the largest change of command (Lemma "ramp disturbance").
module VelocityCartPendulum

import LazySets
import MathematicalSystems
import StaticArrays: SVector, SMatrix
import Dionysos
const PR = Dionysos.Problem
const OPDS = Dionysos.Optim.DiscreteSystems

Base.@kwdef struct Params{T}
    M::T = 1.0  # cart mass
    m::T = 0.2  # pole mass, at the tip
    l::T = 1.0
    g::T = 9.81
end

function dynamics(P::Params)
    m, l, g = P.m, P.l, P.g
    return function (z, v)
        s, c = sincos(z[1])
        w = v[1]
        return SVector(
            (z[2] - m * l * c * w) / (m * l^2),
            -(s / l) * w * z[2] + m * s * c * w^2 + m * g * l * s,
        )
    end
end

angular_velocity(P::Params, θ, p, v) = (p - P.m * P.l * cos(θ) * v) / (P.m * P.l^2)
force(P::Params, θ, θdot, a) =
    (P.M + P.m * sin(θ)^2) * a + P.m * sin(θ) * (P.g * cos(θ) - P.l * θdot^2)

"""
    jacobian_bound(P; p_max)

Bound on the Jacobian of the pole dynamics in (θ, p) for the constant input `u`, valid for
`|p| ≤ p_max`: the diagonal bounds the signed derivative, the off-diagonal its absolute value.
"""
function jacobian_bound(P::Params; p_max)
    m, l, g = P.m, P.l, P.g
    return function (u)
        w = abs(u[1])
        return SMatrix{2, 2}(
            w / l,
            w * p_max / l + m * w^2 + m * g * l,
            1 / (m * l^2),
            w / l,
        )
    end
end

"""
    ramp_disturbance(P; p_max, v_max, Δ_max)

Half-width `b Δ_max` of the additive disturbance that covers a ramp of the input by at most
`Δ_max`, with `b = sup |∂f/∂w|` over `|p| ≤ p_max`, `|w| ≤ v_max`.
"""
ramp_disturbance(P::Params; p_max, v_max, Δ_max) =
    SVector(1 / P.l, p_max / P.l + P.m * v_max) .* Δ_max

"""
    system(P; p_max, v_max, Δ_max)

The pole dynamics in (θ, p), θ ∈ [−π, π] periodic, with the ramp disturbance.
"""
function system(P::Params; p_max, v_max, Δ_max)
    f = dynamics(P)
    z̄ = ramp_disturbance(P; p_max, v_max, Δ_max)
    return MathematicalSystems.NoisyConstrainedBlackBoxControlContinuousSystem(
        (z, u, w) -> f(z, u) + w,
        2,
        1,
        2,
        LazySets.Hyperrectangle(; low = SVector(-π, -p_max), high = SVector(π, p_max)),
        LazySets.Hyperrectangle(; low = SVector(-v_max), high = SVector(v_max)),
        LazySets.Hyperrectangle(; low = -z̄, high = z̄),
    )
end

"""
    nominal(P; p_max, v_max)

The pole dynamics without the ramp disturbance, i.e. under velocity steps.
"""
function nominal(P::Params; p_max, v_max)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamics(P),
        2,
        1,
        LazySets.Hyperrectangle(; low = SVector(-π, -p_max), high = SVector(π, p_max)),
        LazySets.Hyperrectangle(; low = SVector(-v_max), high = SVector(v_max)),
    )
end

# ------------------------------------------------------------------ interval arithmetic
# Just enough to bound the force over a box of states: intervals are (lo, hi) tuples.
_add(a, b) = (a[1] + b[1], a[2] + b[2])
_scale(k, a) = k >= 0 ? (k * a[1], k * a[2]) : (k * a[2], k * a[1])
_mul(a, b) = extrema((a[1] * b[1], a[1] * b[2], a[2] * b[1], a[2] * b[2]))
_sq(a) = a[1] <= 0 <= a[2] ? (0.0, max(a[1]^2, a[2]^2)) : extrema((a[1]^2, a[2]^2))
_within(a, b) = b[1] <= a[1] && a[2] <= b[2]

function _sin(a)
    a[2] - a[1] >= 2π && return (-1.0, 1.0)
    vals = [sin(a[1]), sin(a[2])]
    k = ceil((a[1] - π / 2) / π)
    while π / 2 + k * π <= a[2]
        push!(vals, sin(π / 2 + k * π))
        k += 1
    end
    return extrema(vals)
end
_cos(a) = _sin((a[1] + π / 2, a[2] + π / 2))

function _rates(P::Params, Θ, Pp, W)
    m, l, g = P.m, P.l, P.g
    S, C = _sin(Θ), _cos(Θ)
    θdot = _scale(1 / (m * l^2), _add(Pp, _scale(-m * l, _mul(C, W))))
    pdot = _add(
        _add(_scale(-1 / l, _mul(_mul(S, W), Pp)), _scale(m, _mul(_mul(S, C), _sq(W)))),
        _scale(m * g * l, S),
    )
    return θdot, pdot
end

"""
    tube(P, Θ0, P0, W, τ)

A box enclosing the solutions from the box `Θ0 × P0` over `[0, τ]` under any input in the
interval `W`: an interval Picard iteration, accepted once the box maps into itself.
"""
function tube(P::Params, Θ0, P0, W, τ)
    B = ((Θ0[1] - 0.5, Θ0[2] + 0.5), (P0[1] - 0.5, P0[2] + 0.5))
    found = nothing
    for _ in 1:30
        θdot, pdot = _rates(P, B[1], B[2], W)
        next = (
            (Θ0[1] + τ * min(0.0, θdot[1]), Θ0[2] + τ * max(0.0, θdot[2])),
            (P0[1] + τ * min(0.0, pdot[1]), P0[2] + τ * max(0.0, pdot[2])),
        )
        if _within(next[1], B[1]) && _within(next[2], B[2])
            found = next
            B = next
        elseif found === nothing
            B = (
                (min(B[1][1], next[1][1]) - 0.5, max(B[1][2], next[1][2]) + 0.5),
                (min(B[2][1], next[2][1]) - 0.5, max(B[2][2], next[2][2]) + 0.5),
            )
        else
            break
        end
    end
    found === nothing && error("no a-priori enclosure found")
    return found
end

"""
    admissible(P, τ, F_max, Δ_max, Θ0, P0, v, u)

The torque admissibility of the paper for ramps over the whole step: the command `u` may follow
`v` from the box `Θ0 × P0` if `|u − v| ≤ Δ_max` and the force
`(M + m sin² θ)(u − v)/τ + m sin θ (g cos θ − l θ̇²)` stays within `±F_max` over the tube of the step.
"""
function admissible(P::Params, τ, F_max, Δ_max, Θ0, P0, v, u)
    abs(u - v) <= Δ_max + 1e-9 || return false
    W = (min(v, u), max(v, u))
    Θ, Pp = tube(P, Θ0, P0, W, τ)
    θdot, _ = _rates(P, Θ, Pp, W)
    S, C = _sin(Θ), _cos(Θ)
    α = _add((P.M, P.M), _scale(P.m, _sq(S)))
    β = _scale(P.m, _mul(S, _add(_scale(P.g, C), _scale(-P.l, _sq(θdot)))))
    F = _add(_scale((u - v) / τ, α), β)
    return -F_max <= F[1] && F[2] <= F_max
end

"""
    ramp_sensitivity(P, Θ, Pp, W)

Componentwise bound on `|∂f/∂w|` over the box `Θ × Pp` and the inputs `W`: the cell-wise `b` of
the ramp disturbance.
"""
function ramp_sensitivity(P::Params, Θ, Pp, W)
    S, C = _sin(Θ), _cos(Θ)
    d1 = _scale(1 / P.l, C)
    d2 = _add(_scale(-1 / P.l, _mul(S, Pp)), _scale(P.m, _mul(_mul(S, C), _scale(2.0, W))))
    return SVector(max(abs(d1[1]), abs(d1[2])), max(abs(d2[1]), abs(d2[2])))
end

"""
    torque_limit(P, τ, F_max, Δ_max; target_input = 0)

[`admissible`](@ref) as a constraint on consecutive inputs, written on concrete cells, for the
grid abstraction. Starts at rest; ends at rest unless `target_input = nothing`.
"""
function torque_limit(P::Params, τ, F_max, Δ_max; target_input = SVector(0.0))
    check = function (cell, u_prev, u)
        lo, hi = LazySets.low(cell), LazySets.high(cell)
        return admissible(
            P,
            τ,
            F_max,
            Δ_max,
            (lo[1], hi[1]),
            (lo[2], hi[2]),
            u_prev[1],
            u[1],
        )
    end
    return OPDS.StateInputVariation(check; target_input, initial_input = SVector(0.0))
end

"""
    swing_up_problem(system; θ0, θ_target, p_target)

From rest near the bottom (θ = π) to a neighbourhood of the upright, `|θ| ≤ θ_target`,
`|p| ≤ p_target`, from which the cart can be brought to rest.
"""
function swing_up_problem(system; θ0, θ_target, p_target)
    _I_ = LazySets.Hyperrectangle(;
        low = SVector(θ0 - 1e-3, -1e-3),
        high = SVector(θ0 + 1e-3, 1e-3),
    )
    _T_ = LazySets.Hyperrectangle(;
        low = SVector(-θ_target, -p_target),
        high = SVector(θ_target, p_target),
    )
    return PR.OptimalControlProblem(system, _I_, _T_, nothing, (x, u) -> 1.0, PR.Infinity())
end

end # module
