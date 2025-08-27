
# """
#     TimeSymbolicModel{N, T}

# Structure representing the symbolic abstraction of time for a mode in a hybrid system.

# # Fields
# - `tsteps::SVector{N, Float64}`: Discrete time steps
# - `time_domain::T`: The time domain (e.g., a HyperRectangle)  
# - `is_time_active::Bool`: True if time evolves, false if time is frozen

# # Type Parameters
# - `N`: Number of time steps 
# - `T`: Type of the time domain
# """
struct TimeSymbolicModel{N, T}
    tsteps::SVector{N, Float64}
    time_domain::T
    is_time_active::Bool
end

# """
#     TimeSymbolicModel(sys::MathematicalSystems.ConstrainedLinearContinuousSystem, tstep::Float64)

# Construct a symbolic time model from a continuous system and a time discretization step.

# # Time Evolution Logic
# - If A ≈ I (identity): time evolves, generates grid from tmin to tmax with tstep
# - If A ≈ 0 (zero): time is frozen, single time step at 0.0

# # Arguments
# - `sys`: 1D ConstrainedLinearContinuousSystem representing time dynamics
# - `tstep::Float64`: Time discretization step (ignored if time is frozen)

# # Returns
# - `TimeSymbolicModel`: Parameterized by number of steps and domain type
# """
function TimeSymbolicModel(
    sys::MathematicalSystems.ConstrainedLinearContinuousSystem,
    tstep::Float64,
)
    A = sys.A
    X = sys.X

    # Extract time domain bounds
    tmin, tmax = try
        X.lb[1], X.ub[1]
    catch
        error("Time domain X must have .lb and .ub fields (e.g., UT.HyperRectangle)")
    end

    # Efficient matrix type detection
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

# Efficient matrix type checking helpers
_is_identity_matrix(A::AbstractMatrix)::Bool = A ≈ I(size(A, 1))
_is_zero_matrix(A::AbstractMatrix)::Bool = all(iszero, A)

# """
#     floor_time2int(tm::TimeSymbolicModel, t::Float64) -> Int

# Return the index of the largest time step ≤ t.
# If time is frozen, always returns 1.
# """
function floor_time2int(tm::TimeSymbolicModel, t::Float64)::Int
    if tm.is_time_active
        idx = findlast(x -> x <= t, tm.tsteps)
        return idx === nothing ? 1 : idx
    else
        return 1
    end
end

# Alias for backward compatibility
const time2int = floor_time2int

# """
#     int2time(tm::TimeSymbolicModel, idx::Int) -> Float64

# Return the time value corresponding to index idx.
# If time is frozen, always returns 0.0.
# """
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
