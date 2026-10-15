# Velocity-controlled pendulum: a velocity servo drives the joint, so the commanded angular velocity
# u is the input, held over each sampling period (step commands), and the model is the
# one-dimensional integrator `θ̇ = u` on the circle.
#
# Gravity has left the dynamics: it reappears only in the torque the servo must produce, in units of
# m l² as in `SimplePendulum` (θ from the bottom). Between steps the servo holds the velocity
# against gravity with the torque (g/l) sin θ, and a step from v to u needs the impulse u − v.
module VelocityPendulum

import LazySets
import MathematicalSystems
import StaticArrays: SVector
import Dionysos
const PR = Dionysos.Problem
const MP = Dionysos.Mapping

"""
    discretization(; τ = 0.1, n = 125, levels = 4)

Commensurable grids on the circle: `n` cells of width `h = 2π/n`, and speed step `δu = h/τ`, so
that `θ̇ = u` maps cells exactly onto cells, including across the wrap at ±π. The cells tile
`[−π, π)` from −π; with `n` odd, θ = 0 is a cell center. The commands are the `2 levels + 1`
multiples of `δu` up to `u_max = levels δu`.
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

The velocity-input pendulum `θ̇ = u` on the full circle, with `|u| ≤ disc.u_max`; abstract it with
a periodic mapping of period 2π starting at −π.
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
    swing_up_problem(system; θ0 = 0.0, target_halfwidth = 15π/180, h = 0.05)

From `θ0` (a set of width `h/2` around it, inside one cell of width `h`) to the top,
`|θ − π| ≤ target_halfwidth`, at unit cost per step.
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
