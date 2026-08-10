"""
    ClockAbstraction{N, T}

Symbolic abstraction of a 1-D monotone clock (the time axis of one mode of a timed
hybrid system, or a standalone clock lifted onto a continuous model).

`tsteps` are the discretized time values; `domain` is the concrete time domain `X`;
`is_active` distinguishes an evolving clock (`ẋ = 1`, matrix `A = I`) from a frozen
one (`A = 0`).
"""
struct ClockAbstraction{N, T}
    tsteps::SVector{N, Float64}
    domain::T
    is_active::Bool
end

"""
    ClockAbstraction(sys::ConstrainedLinearContinuousSystem, tstep)

Build a [`ClockAbstraction`](@ref) from a linear time subsystem. The dynamics
matrix `A` must be the identity (time evolves, discretized with step `tstep` over
the box time domain) or zero (time frozen, a single step at `0.0`).
"""
function ClockAbstraction(sys::MS.ConstrainedLinearContinuousSystem, tstep::Float64)
    A = sys.A
    X = sys.X

    tmin, tmax = try
        LazySets.low(X, 1), LazySets.high(X, 1)
    catch
        error("Time domain X must have box bounds (e.g., a `LazySets.Hyperrectangle`)")
    end

    if _is_identity_matrix(A)  # Time evolves
        tsteps_vec = collect(tmin:tstep:tmax)
        N = length(tsteps_vec)
        tsteps = SVector{N, Float64}(tsteps_vec)
        return ClockAbstraction{N, typeof(X)}(tsteps, X, true)
    elseif _is_zero_matrix(A)  # Time is frozen
        tsteps = SVector{1, Float64}(0.0)
        return ClockAbstraction{1, typeof(X)}(tsteps, X, false)
    else
        error("Matrix A must be identity (time active) or zero (time frozen). Got: $A")
    end
end

"""
    ClockAbstraction(domain, tstep::Float64)

Build an *active* [`ClockAbstraction`](@ref) directly from a 1-D time `domain`
(a box interval) discretized with step `tstep`. Used to lift a continuous model
with a clock (no time subsystem required).
"""
function ClockAbstraction(domain, tstep::Float64)
    tmin, tmax = LazySets.low(domain, 1), LazySets.high(domain, 1)
    tsteps_vec = collect(tmin:tstep:tmax)
    N = length(tsteps_vec)
    return ClockAbstraction{N, typeof(domain)}(
        SVector{N, Float64}(tsteps_vec),
        domain,
        true,
    )
end

_is_identity_matrix(A::AbstractMatrix)::Bool = A ≈ LA.I(size(A, 1))
_is_zero_matrix(A::AbstractMatrix)::Bool = all(iszero, A)

"""
    floor_time2int(tm::ClockAbstraction, t) -> Int

Index of the largest time step `≤ t` (always `1` when time is frozen).
"""
function floor_time2int(tm::ClockAbstraction, t::Float64)::Int
    if tm.is_active
        idx = findlast(x -> x <= t, tm.tsteps)
        return idx === nothing ? 1 : idx
    else
        return 1
    end
end

"""
    ceil_time2int(tm::ClockAbstraction, t) -> Int

Index of the smallest time step `≥ t` (always `1` when time is frozen).
"""
function ceil_time2int(tm::ClockAbstraction, t::Float64)::Int
    if tm.is_active
        idx = findfirst(x -> x >= t, tm.tsteps)
        return idx === nothing ? length(tm.tsteps) : idx
    else
        return 1
    end
end

"""
    int2time(tm::ClockAbstraction, idx) -> Float64

Time value at index `idx` (always `0.0` when time is frozen).
"""
function int2time(tm::ClockAbstraction, idx::Int)::Float64
    if tm.is_active
        @inbounds return tm.tsteps[idx]
    else
        return 0.0
    end
end
