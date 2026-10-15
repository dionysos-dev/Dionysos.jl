# Velocity-controlled cart-pole: a velocity servo drives the cart, whose commanded velocity v is the
# input, held over each sampling period (step commands). With the pole angle θ, measured from the
# upright, and the momentum p = m l² θ̇ + m l cos θ v conjugate to θ, a point-mass pole without
# friction obeys
#   θ̇ = (p − m l cos θ v) / (m l²)
#   ṗ = −(sin θ / l) v p + m sin θ cos θ v² + m g l sin θ.
# The cart position enters neither equation, so (θ, p) is the whole state on an unbounded track.
# The cart force that realizes the motion with cart acceleration a = v̇ is
#   F = (M + m sin² θ) a + m sin θ (g cos θ − l θ̇²):
# between steps a = 0 and F is the hold force, and a step from v to u is the impulse
# (M + m sin² θ)(u − v).
module VelocityCartPendulum

import LazySets
import MathematicalSystems
import StaticArrays: SVector, SMatrix
import Dionysos
const PR = Dionysos.Problem

Base.@kwdef struct Params{T}
    M::T = 1.0  # cart mass
    m::T = 0.2  # pole mass, at the tip
    l::T = 1.0
    g::T = 9.81
end

"""
    dynamics(P) -> f

The pole dynamics `f(z, v)` in `z = (θ, p)` under the cart velocity `v`.
"""
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

"""
    angular_velocity(P, θ, p, v)

The angular velocity θ̇ of the pole at the state `(θ, p)` under the cart velocity `v`.
"""
angular_velocity(P::Params, θ, p, v) = (p - P.m * P.l * cos(θ) * v) / (P.m * P.l^2)

"""
    force(P, θ, θdot, a)

The cart force that realizes the cart acceleration `a` at the angle `θ` and angular velocity
`θdot`; with `a = 0`, the hold force.
"""
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
    system(P; p_max, v_max)

The pole dynamics in (θ, p), θ ∈ [−π, π] periodic and `|p| ≤ p_max`, under cart velocities
`|v| ≤ v_max`.
"""
function system(P::Params; p_max, v_max)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamics(P),
        2,
        1,
        LazySets.Hyperrectangle(; low = SVector(-π, -p_max), high = SVector(π, p_max)),
        LazySets.Hyperrectangle(; low = SVector(-v_max), high = SVector(v_max)),
    )
end

"""
    swing_up_problem(system; θ0, θ_target, p_target)

From rest at `θ0`, the pole hanging down for `θ0 = π`, to a neighborhood of the upright,
`|θ| ≤ θ_target` and `|p| ≤ p_target`, at unit cost per step.
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

# ------------------------------------------------------------------------ interval arithmetic
"""
Interval arithmetic, just enough to bound the dynamics and the cart force over a box of states.
An interval is a tuple `(lo, hi)`:

- `iadd(a, b)`, `iscale(k, a)`, `imul(a, b)`, `isq(a)`: sum, scaling by a number, product, square;
- `isin(a)`, `icos(a)`: the exact range of `sin` and `cos` over `a`;
- `iwithin(a, b)`: whether `a ⊆ b`.
"""
iadd(a, b) = (a[1] + b[1], a[2] + b[2])
iscale(k, a) = k >= 0 ? (k * a[1], k * a[2]) : (k * a[2], k * a[1])
imul(a, b) = extrema((a[1] * b[1], a[1] * b[2], a[2] * b[1], a[2] * b[2]))
isq(a) = a[1] <= 0 <= a[2] ? (0.0, max(a[1]^2, a[2]^2)) : extrema((a[1]^2, a[2]^2))
iwithin(a, b) = b[1] <= a[1] && a[2] <= b[2]
function isin(a)
    a[2] - a[1] >= 2π && return (-1.0, 1.0)
    vals = [sin(a[1]), sin(a[2])]
    k = ceil((a[1] - π / 2) / π)
    while π / 2 + k * π <= a[2]
        push!(vals, sin(π / 2 + k * π))
        k += 1
    end
    return extrema(vals)
end
icos(a) = isin((a[1] + π / 2, a[2] + π / 2))

end # module
