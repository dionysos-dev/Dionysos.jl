import LinearAlgebra as LA
using StaticArrays

"""
TimeGridMapping:
- 1D mapping over time values t
- state id q ∈ 1..K
- coord by state is the discrete time value tm.tsteps[q]
- coord to state can be floor/ceil depending on how you want to abstract
"""
struct TimeGridMapping{Tdom} <: AbstractMapping{1,Float64}
    tsteps::Vector{Float64}
    time_domain::Tdom
    is_time_active::Bool
end

# ---- constructors ----

_is_identity_matrix(A::AbstractMatrix)::Bool = A ≈ LA.I(size(A, 1))
_is_zero_matrix(A::AbstractMatrix)::Bool = all(iszero, A)

function TimeGridMapping(
    sys::MathematicalSystems.ConstrainedLinearContinuousSystem,
    tstep::Float64,
)
    A = sys.A
    X = sys.X

    # expect a rectangle-like with lb/ub
    tmin, tmax = try
        X.lb[1], X.ub[1]
    catch
        error("Time domain X must have .lb and .ub (e.g. UT.HyperRectangle)")
    end

    if _is_identity_matrix(A)  # time evolves
        tsteps = collect(tmin:tstep:tmax)
        isempty(tsteps) && error("Empty time grid: check tmin/tmax/tstep.")
        return TimeGridMapping(tsteps, X, true)

    elseif _is_zero_matrix(A)  # discrete time
        return TimeGridMapping([0.0], X, false)

    else
        error("Matrix A must be identity (time active) or zero (time frozen). Got: $A")
    end
end

# ---- required mapping API (minimal) ----

get_n_state(m::TimeGridMapping) = length(m.tsteps)
enum_states(m::TimeGridMapping) = 1:get_n_state(m)

get_coord_by_state(m::TimeGridMapping, q::Int) = begin
    (1 <= q <= get_n_state(m)) || throw(DomainError(q, "time state out of range"))
    return SVector{1,Float64}(m.tsteps[q])
end

"""
Default: floor abstraction (largest tgrid <= t).
If time frozen => always return 1.
"""
function get_state_by_coord(m::TimeGridMapping, t)
    tt = t isa Number ? Float64(t) : Float64(t[1])
    return floor_time2int(m, tt)
end

# ---- Additionnal Helpers ----

function floor_time2int(m::TimeGridMapping, t::Float64)::Int
    if m.is_time_active
        idx = findlast(x -> x <= t, m.tsteps)
        return idx === nothing ? 1 : idx
    else
        return 1
    end
end

function ceil_time2int(m::TimeGridMapping, t::Float64)::Int
    if m.is_time_active
        idx = findfirst(x -> x >= t, m.tsteps)
        return idx === nothing ? length(m.tsteps) : idx
    else
        return 1
    end
end

function int2time(m::TimeGridMapping, q::Int)::Float64
    if m.is_time_active
        (1 <= q <= get_n_state(m)) || throw(DomainError(q, "time state out of range"))
        return m.tsteps[q]
    else
        return 0.0
    end
end