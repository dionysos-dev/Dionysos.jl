
"""
    TimeSymbolicModel

Structure representing the symbolic abstraction of time for a mode in a hybrid system.

# Fields
- `tsteps::Vector{Float64}`: Discrete time steps
- `time_domain::Any`: The time domain (e.g., a HyperRectangle)
- `is_active::Bool`: True if time evolves, false if time is frozen
"""
struct TimeSymbolicModel
    tsteps::Vector{Float64}
    time_domain::Any
    is_active::Bool
end

"""
    TimeSymbolicModel(sys::MathematicalSystems.ConstrainedLinearContinuousSystem, tstep::Float64)

Construct a symbolic time model from a continuous system and a time discretization step.
If the system is time-evolving (A = I), generates a grid of time steps.
If the system is time-frozen (A = 0), only one time step (0.0) is used.

# Arguments
- `sys`: The time system (should be a 1D ConstrainedLinearContinuousSystem)
- `tstep::Float64`: Time discretization step

# Returns
- `TimeSymbolicModel`: The symbolic time model
"""
function BuildTimeSymbolicModel(
    sys::MathematicalSystems.ConstrainedLinearContinuousSystem,
    tstep::Float64,
)
    A = sys.A
    X = sys.X
    # Extract time domain bounds (assumes X has .lb and .ub fields)
    tmin, tmax = try
        X.lb[1], X.ub[1]
    catch
        error(
            "Time domain X must have .lb and .ub fields (e.g., UT.HyperRectangle from Dionysos).",
        )
    end
    if A == ones(size(A))  # Identity matrix: time evolves
        tsteps = collect(tmin:tstep:tmax)
        return TimeSymbolicModel(tsteps, X, true)
    elseif A == zeros(size(A))  # Zero matrix: time is frozen
        tsteps = [0.0]
        return TimeSymbolicModel(tsteps, X, false)
    else
        error("Matrix A must be 0 or 1 for time handling.")
    end
end

"""
    time2int(tm::TimeSymbolicModel, t::Real)

Return the index of the largest time step less than or equal to t.
If time is frozen, always returns 1.
"""
function time2int(tm::TimeSymbolicModel, t::Real)
    if tm.is_active
        idx = findlast(x -> x <= t, tm.tsteps)
        return idx === nothing ? 1 : idx
    else
        return 1
    end
end

"""
    int2time(tm::TimeSymbolicModel, idx::Int)

Return the time value corresponding to a given index.
If time is frozen, always returns 0.0.
"""
function int2time(tm::TimeSymbolicModel, idx::Int)
    if tm.is_active
        return tm.tsteps[idx]
    else
        return 0.0
    end
end

"""
    ceil_time2int(tm::TimeSymbolicModel, t::Real)

Return the index of the smallest time step greater than or equal to t.
If time is frozen, always returns 1.
"""
function ceil_time2int(tm::TimeSymbolicModel, t::Real)
    if tm.is_active
        idx = findfirst(x -> x >= t, tm.tsteps)
        return idx === nothing ? length(tm.tsteps) : idx
    else
        return 1
    end
end

"""
    floor_time2int(tm::TimeSymbolicModel, t::Real)

Return the index of the largest time step less than or equal to t.
If time is frozen, always returns 1.
"""
function floor_time2int(tm::TimeSymbolicModel, t::Real)
    if tm.is_active
        idx = findlast(x -> x <= t, tm.tsteps)
        return idx === nothing ? 1 : idx
    else
        return 1
    end
end
