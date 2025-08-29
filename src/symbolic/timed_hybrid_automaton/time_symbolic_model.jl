# Structure representing the symbolic abstraction of time for a mode in a hybrid system.
struct TimeSymbolicModel{N, T}
    tsteps::SVector{N, Float64}
    time_domain::T
    is_time_active::Bool
end

function TimeSymbolicModel(
    sys::MathematicalSystems.ConstrainedLinearContinuousSystem,
    tstep::Float64,
)
    A = sys.A
    X = sys.X

    tmin, tmax = try
        X.lb[1], X.ub[1]
    catch
        error("Time domain X must have .lb and .ub fields (e.g., UT.HyperRectangle)")
    end

    if _is_identity_matrix(A)  # Time evolves
        tsteps_vec = collect(tmin:tstep:tmax)
        N = length(tsteps_vec)
        tsteps = SVector{N, Float64}(tsteps_vec)
        return TimeSymbolicModel{N, typeof(X)}(tsteps, X, true)
    elseif _is_zero_matrix(A)  # Time is frozen
        tsteps = SVector{1, Float64}(0.0)
        return TimeSymbolicModel{1, typeof(X)}(tsteps, X, false)
    else
        error("Matrix A must be identity (time active) or zero (time frozen). Got: $A")
    end
end

_is_identity_matrix(A::AbstractMatrix)::Bool = A ≈ I(size(A, 1))
_is_zero_matrix(A::AbstractMatrix)::Bool = all(iszero, A)

# floor_time2int(tm::TimeSymbolicModel, t::Float64) -> Int
# Return the index of the largest time step ≤ t.
# If time is frozen, always returns 1.
function floor_time2int(tm::TimeSymbolicModel, t::Float64)::Int
    if tm.is_time_active
        idx = findlast(x -> x <= t, tm.tsteps)
        return idx === nothing ? 1 : idx
    else
        return 1
    end
end


# int2time(tm::TimeSymbolicModel, idx::Int) -> Float64
# Return the time value corresponding to index idx.
# If time is frozen, always returns 0.0.
function int2time(tm::TimeSymbolicModel, idx::Int)::Float64
    if tm.is_time_active
        @inbounds return tm.tsteps[idx]
    else
        return 0.0
    end
end

# """
#     ceil_time2int(tm::TimeSymbolicModel, t::Float64) -> Int

# Return the index of the smallest time step ≥ t.
# If time is frozen, always returns 1.
# """
function ceil_time2int(tm::TimeSymbolicModel, t::Float64)::Int
    if tm.is_time_active
        idx = findfirst(x -> x >= t, tm.tsteps)
        return idx === nothing ? length(tm.tsteps) : idx
    else
        return 1
    end
end
