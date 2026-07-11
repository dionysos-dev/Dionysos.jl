# ----------------------------
# Time horizon
# ----------------------------

"""
    Infinity

Sentinel for an infinite time horizon. A problem whose `time` is `Infinity()`
must be satisfied over an unbounded horizon (e.g. an infinite-horizon safety or
reach-and-stay specification). It is never finite, so solvers branch on
`isfinite(problem.time)` to separate the finite from the infinite-horizon case.
"""
struct Infinity end

Base.isfinite(::Infinity) = false
Base.isinf(::Infinity) = true

"""
    discretize_time(time, Δt::Real; round_up = true)

Convert a continuous-time horizon `time` into a number of discrete steps of
duration `Δt`. `round_up` selects the conservative direction: `true` rounds up
(used by "for at least `T`" specifications such as safety and reach-and-stay),
`false` rounds down (used by "within at most `T`" specifications such as
reach-avoid). An `Infinity()` horizon stays `Infinity()`.
"""
discretize_time(time::Real, Δt::Real; round_up = true) =
    round_up ? ceil(Int, time / Δt) : floor(Int, time / Δt)
discretize_time(::Infinity, Δt::Real; round_up = true) = Infinity()
