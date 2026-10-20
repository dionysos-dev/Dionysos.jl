# Frictionless cart-pole with a point-mass pole at the tip, θ measured from the upright, in two
# models. The cart position enters neither the dynamics of the pole nor the tasks on an unbounded
# track, so both leave it out.
#
# - Velocity model (F): a velocity servo drives the cart, whose commanded velocity v is held over
#   each sampling period. With the momentum p = m l² θ̇ + m l cos θ v conjugate to θ, the state is
#   (θ, p) and
#     θ̇ = (p − m l cos θ v) / (m l²)
#     ṗ = −(sin θ / l) v p + m sin θ cos θ v² + m g l sin θ.
#   The cart force that realizes the motion with cart acceleration a = v̇ is
#     F = (M + m sin² θ) a + m sin θ (g cos θ − l θ̇²):
#   between steps a = 0 and F is the hold force, and a step from v to u is the impulse
#   (M + m sin² θ)(u − v).
# - Force model (T): the cart force F is the input, held over each sampling period, and the state is
#   (θ, ẋ, θ̇). Its growth bound uses a Jacobian bound per cell, over an a-priori enclosure of the
#   step from the cell.
# - Velocity servo: a PI loop on the force model that tracks a commanded cart velocity under a
#   force limit, to replay a controller of the velocity model without its ideal servo.
#
# The bounds use the interval arithmetic below, on (lo, hi) tuples; like the integration, they are
# sound up to floating-point rounding.
module PointMassCartPendulum

import LazySets
import MathematicalSystems
import StaticArrays: SVector, SMatrix
import Dionysos
const PR = Dionysos.Problem
const ST = Dionysos.System
const OPDS = Dionysos.Optim.DiscreteSystems

Base.@kwdef struct Params{T}
    M::T = 1.0  # cart mass
    m::T = 0.2  # pole mass, at the tip
    l::T = 1.0
    g::T = 9.81
end

# ------------------------------------------------------------------------ interval arithmetic
"""
Interval arithmetic, just enough to bound the dynamics and the cart force over a box of states.
An interval is a tuple `(lo, hi)`:

- `iadd(a, b)`, `iscale(k, a)`, `imul(a, b)`, `isq(a)`, `idiv(a, d)`: sum, scaling by a number,
  product, square, and quotient by an interval `d > 0`;
- `isin(a)`, `icos(a)`: the exact range of `sin` and `cos` over `a`;
- `ihull0(a)`: the hull of `a` and 0; `iabsmax(a)`: the largest absolute value in `a`;
- `iwithin(a, b)`: whether `a ⊆ b`.
"""
iadd(a, b) = (a[1] + b[1], a[2] + b[2])
iscale(k, a) = k >= 0 ? (k * a[1], k * a[2]) : (k * a[2], k * a[1])
imul(a, b) = extrema((a[1] * b[1], a[1] * b[2], a[2] * b[1], a[2] * b[2]))
isq(a) = a[1] <= 0 <= a[2] ? (0.0, max(a[1]^2, a[2]^2)) : extrema((a[1]^2, a[2]^2))
idiv(a, d) = extrema((a[1] / d[1], a[1] / d[2], a[2] / d[1], a[2] / d[2]))
ihull0(a) = (min(0.0, a[1]), max(0.0, a[2]))
iabsmax(a) = max(abs(a[1]), abs(a[2]))
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

# The interval of coordinate `i` of a box.
interval(box, i) = (LazySets.low(box, i), LazySets.high(box, i))

# ------------------------------------------------------------------------- velocity model (F)
"""
    velocity_dynamics(P) -> f

The pole dynamics `f(z, v)` in `z = (θ, p)` under the cart velocity `v`.
"""
function velocity_dynamics(P::Params)
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
    cart_force(P, θ, θdot, a)

The cart force that realizes the cart acceleration `a` at the angle `θ` and angular velocity
`θdot`; with `a = 0`, the hold force.
"""
cart_force(P::Params, θ, θdot, a) =
    (P.M + P.m * sin(θ)^2) * a + P.m * sin(θ) * (P.g * cos(θ) - P.l * θdot^2)

"""
    velocity_jacobian_bound(P; p_max) -> u -> SMatrix

Bound on the Jacobian of the pole dynamics in (θ, p) for the constant input `u`, valid for
`|p| ≤ p_max`: the diagonal bounds the signed derivative, the off-diagonal its absolute value.
"""
function velocity_jacobian_bound(P::Params; p_max)
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
    velocity_system(P; p_max, v_max)

The pole dynamics in (θ, p), θ ∈ [−π, π] periodic and `|p| ≤ p_max`, under cart velocities
`|v| ≤ v_max`.
"""
function velocity_system(P::Params; p_max, v_max)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        velocity_dynamics(P),
        2,
        1,
        LazySets.Hyperrectangle(; low = SVector(-π, -p_max), high = SVector(π, p_max)),
        LazySets.Hyperrectangle(; low = SVector(-v_max), high = SVector(v_max)),
    )
end

# Upper bound of the inertia α(θ) = M + m sin²θ over an interval of angles.
inertia_bound(P::Params, Θ) = P.M + P.m * isq(isin(Θ))[2]

"""
    impulse_budget(P, I_max)

The impulse budget |α(θ)(u − v)| ≤ I_max as a constraint on consecutive cart velocities, carried
to a cell by the upper bound of α over its angles. The motion starts at rest; any velocity is
allowed on arrival.
"""
function impulse_budget(P::Params, I_max)
    compatible(cell, v, u) =
        inertia_bound(P, interval(cell, 1)) * abs(u[1] - v[1]) <= I_max + 1e-9
    return OPDS.StateInputVariation(
        compatible;
        target_input = nothing,
        initial_input = SVector(0.0),
    )
end

# ---------------------------------------------------------------------------- force model (T)
"""
    force_dynamics(P) -> f

The dynamics `f(y, F)` in `y = (θ, ẋ, θ̇)` under the cart force `F`.
"""
function force_dynamics(P::Params)
    return function (y, u)
        F = u[1]
        s, c = sincos(y[1])
        a = (F + P.m * s * (P.l * y[3]^2 - P.g * c)) / (P.M + P.m * s^2)
        return SVector(y[3], a, (P.g * s - c * a) / P.l)
    end
end

"""
    force_system(P; v_max, ω_max, F_max)

The force model in (θ, ẋ, θ̇), θ ∈ [−π, π] periodic, `|ẋ| ≤ v_max`, `|θ̇| ≤ ω_max`, under cart
forces `|F| ≤ F_max`.
"""
function force_system(P::Params; v_max, ω_max, F_max)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        force_dynamics(P),
        3,
        1,
        LazySets.Hyperrectangle(;
            low = SVector(-π, -v_max, -ω_max),
            high = SVector(π, v_max, ω_max),
        ),
        LazySets.Hyperrectangle(; low = SVector(-F_max), high = SVector(F_max)),
    )
end

# Interval rates of (θ, ẋ, θ̇) over the box Θ × Ω (ẋ does not enter), and the pieces the Jacobian
# reuses.
function force_rates(P::Params, Θ, Ω, F)
    S, C = isin(Θ), icos(Θ)
    D = iadd((P.M, P.M), iscale(P.m, isq(S)))
    N = iadd((F, F), iscale(P.m, imul(S, iadd(iscale(P.l, isq(Ω)), iscale(-P.g, C)))))
    A = idiv(N, D)
    Θdd = iscale(1 / P.l, iadd(iscale(P.g, S), iscale(-1.0, imul(C, A))))
    return A, Θdd, (S, C, D, N)
end

"""
    force_enclosure(P, Y0, F, τ) -> Union{Nothing, NTuple{3}}

A box of (θ, ẋ, θ̇), as three intervals, that contains every solution from the box `Y0` under the
constant force `F` over [0, τ]: an interval Picard iteration `B ⊇ Y0 + [0, τ] f(B)`, started from
`Y0` widened and accepted once the box maps into itself. `nothing` if none is found in 60
iterations.
"""
function force_enclosure(P::Params, Y0, F, τ)
    widen(B, k) = (
        (B[1][1] - k, B[1][2] + k),
        (B[2][1] - k, B[2][2] + k),
        (B[3][1] - 4k, B[3][2] + 4k),
    )
    B = widen(Y0, 0.25)
    found = nothing
    for _ in 1:60
        A, Θdd, _ = force_rates(P, B[1], B[3], F)
        next = (
            iadd(Y0[1], iscale(τ, ihull0(B[3]))),
            iadd(Y0[2], iscale(τ, ihull0(A))),
            iadd(Y0[3], iscale(τ, ihull0(Θdd))),
        )
        if all(iwithin(next[i], B[i]) for i in 1:3)
            found = next
            B = next
        elseif found === nothing
            B = widen(
                (
                    (min(B[1][1], next[1][1]), max(B[1][2], next[1][2])),
                    (min(B[2][1], next[2][1]), max(B[2][2], next[2][2])),
                    (min(B[3][1], next[3][1]), max(B[3][2], next[3][2])),
                ),
                0.25,
            )
        else
            break
        end
    end
    return found
end

"""
    force_jacobian_bound(P, Θ, Ω, F) -> SMatrix

Bound on the Jacobian of the force model over the box `Θ × Ω` of (θ, θ̇) (ẋ does not enter) under
the force `F`: signed upper bounds on the diagonal, absolute values off it.
"""
function force_jacobian_bound(P::Params, Θ, Ω, F)
    _, _, (S, C, D, N) = force_rates(P, Θ, Ω, F)
    m, l, g = P.m, P.l, P.g
    A = idiv(N, D)
    cos2 = icos((2Θ[1], 2Θ[2]))
    num = iadd(
        imul(iscale(m, iadd(iscale(l, imul(C, isq(Ω))), iscale(-g, cos2))), D),
        iscale(-2m, imul(imul(S, C), N)),
    )
    dA_dθ = idiv(num, isq(D))
    dA_dω = idiv(iscale(2m * l, imul(S, Ω)), D)
    dB_dθ =
        iscale(1 / l, iadd(iadd(iscale(g, C), imul(S, A)), iscale(-1.0, imul(C, dA_dθ))))
    dB_dω = iscale(-1 / l, imul(C, dA_dω))
    return SMatrix{3, 3}(
        0.0,
        iabsmax(dA_dθ),
        iabsmax(dB_dθ),
        0.0,
        0.0,
        0.0,
        1.0,
        iabsmax(dA_dω),
        dB_dω[2],
    )
end

"""
    force_cellwise_bound(P) -> Dionysos.System.CellwiseBound

The Jacobian bound of the force model per cell, over the enclosure of the step from the cell
given by [`force_enclosure`](@ref): the growth bound of the force model's abstraction.
"""
function force_cellwise_bound(P::Params)
    enclosure(cell, u, τ) = force_enclosure(
        P,
        (interval(cell, 1), interval(cell, 2), interval(cell, 3)),
        u[1],
        τ,
    )
    bound(E, u) = force_jacobian_bound(P, E[1], E[3], u[1])
    return ST.CellwiseBound(enclosure, bound)
end

"""
    in_force_target(P, Θ, V, Ω; θ_target, p_target)

Whether the box `Θ × V × Ω` of (θ, ẋ, θ̇) lies in the upright target `|θ| ≤ θ_target`,
`|m l² θ̇ + m l cos θ ẋ| ≤ p_target`: the momentum target of (F) in the coordinates of (T).
"""
function in_force_target(P::Params, Θ, V, Ω; θ_target, p_target)
    (-θ_target <= Θ[1] && Θ[2] <= θ_target) || return false
    p = iadd(iscale(P.m * P.l^2, Ω), iscale(P.m * P.l, imul(icos(Θ), V)))
    return -p_target <= p[1] && p[2] <= p_target
end

# ------------------------------------------------------------------------------ velocity servo
"""
    velocity_state(P, y) -> SVector

The state `(θ, p)` of the velocity model at the state `y = (θ, ẋ, θ̇)` of the force model.
"""
velocity_state(P::Params, y) =
    SVector(y[1], P.m * P.l^2 * y[3] + P.m * P.l * cos(y[1]) * y[2])

"""
    force_state(P, z, v) -> SVector

The state `(θ, ẋ, θ̇)` of the force model at the state `z = (θ, p)` of the velocity model, with
the cart velocity `v`.
"""
force_state(P::Params, z, v) = SVector(z[1], v, angular_velocity(P, z[1], z[2], v))

struct VelocityError end
(::VelocityError)(x, _, _) = x[4] - x[2]

"""
    velocity_loop(; Kp, Ki, F_max, dt) -> Dionysos.System.PIDControllers.PIDController

A PI loop that drives the cart velocity of the force model to a setpoint, the servo that the
velocity model idealizes. It runs every `dt`; its input is the state `(θ, ẋ, θ̇)` followed by the
setpoint, and its output the cart force, clamped to `|F| ≤ F_max`, with anti-windup.
"""
function velocity_loop(; Kp, Ki, F_max, dt)
    return ST.PIDControllers.PIDControllerVector(;
        Kp,
        Ki,
        ref = nothing,
        error = VelocityError(),
        dt,
        umin = -F_max,
        umax = F_max,
        e0 = 0.0,
    )
end

end # module
