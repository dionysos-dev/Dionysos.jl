import Dionysos
import StaticArrays: SVector

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const OP = DI.Optim

const DEFAULT_PERIODIC_DIMS = SVector(3, 4)
const DEFAULT_PERIODIC_PERIODS = SVector(2pi, 2pi)
const DEFAULT_PERIODIC_START = SVector(-pi, -pi)

"""
Ensure compatibility with legacy symbolic helpers used by some benchmark scripts.
"""
function ensure_symbolic_format_compat!()
    if !isdefined(DI.Symbolic, :format_input_set)
        @eval DI.Symbolic format_input_set(args...) = $(UT).format_input_set(args...)
    end
    if !isdefined(DI.Symbolic, :format_noise_set)
        @eval DI.Symbolic format_noise_set(args...) = $(UT).format_noise_set(args...)
    end
    return nothing
end

"""
Validate periodic mapping settings for a state space of size nx.
"""
function audit_periodicity!(periodic_dims, periodic_periods, periodic_start; nx::Int)
    nd = length(periodic_dims)
    np = length(periodic_periods)
    ns = length(periodic_start)

    nd == np || error("periodicity invalid: length(periodic_dims)=$nd != length(periodic_periods)=$np")
    nd == ns || error("periodicity invalid: length(periodic_dims)=$nd != length(periodic_start)=$ns")
    nd >= 1 || error("periodicity invalid: expected at least one periodic dimension")

    dims = Int.(collect(periodic_dims))
    allunique(dims) || error("periodicity invalid: periodic_dims contains duplicates: $(dims)")

    for i in eachindex(dims)
        d = dims[i]
        1 <= d <= nx || error("periodicity invalid: dimension $d out of bounds 1:$nx")
        p = Float64(periodic_periods[i])
        p > 0.0 || error("periodicity invalid: periodic_periods[$i]=$p must be > 0")
    end

    return nothing
end

"""
Unwrap periodic state dimensions to remove +/-2pi discontinuities.
"""
function unwrap_periodic_state_list(state_list, periodic_dims, periodic_periods)
    isempty(state_list) && return state_list

    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims and periodic_periods must have the same length")

    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])

        1 <= d <= nx || error("invalid periodic dimension: $d")
        p > 0 || error("invalid period (<= 0): $p")

        for k in 2:length(xs)
            delta_raw = xs[k][d] - xs[k - 1][d]
            delta_min = delta_raw - round(delta_raw / p) * p
            xs[k][d] = xs[k - 1][d] + delta_min
        end
    end

    return [SVector{nx, Float64}(x) for x in xs]
end

"""
Wrap one state vector in-place into periodic ranges.
"""
function wrap_vector_periodic!(x::AbstractVector{Float64}, periodic_dims, periodic_periods, periodic_start)
    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        s = Float64(periodic_start[i])
        x[d] = mod(x[d] - s, p) + s
    end
    return x
end

"""
Wrap all states in a trajectory-like list into periodic ranges.
"""
function wrap_periodic_state_list(state_list, periodic_dims, periodic_periods, periodic_start)
    isempty(state_list) && return state_list

    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims and periodic_periods must have the same length")
    length(periodic_dims) == length(periodic_start) ||
        error("periodic_dims and periodic_start must have the same length")

    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]

    for x in xs
        wrap_vector_periodic!(x, periodic_dims, periodic_periods, periodic_start)
    end

    return [SVector{nx, Float64}(x) for x in xs]
end

"""
Count trajectory jumps larger than half-period on periodic dimensions.
"""
function count_periodic_wrap_jumps(state_list, periodic_dims, periodic_periods)
    length(state_list) <= 1 && return 0

    njump = 0
    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        threshold = 0.5 * p

        for k in 2:length(state_list)
            delta = Float64(state_list[k][d] - state_list[k - 1][d])
            abs(delta) > threshold && (njump += 1)
        end
    end

    return njump
end

"""
Unwrap ellipsoid centers only (matrix P unchanged).
"""
function unwrap_ellipsoid_centers(ellipsoids, periodic_dims, periodic_periods)
    isempty(ellipsoids) && return ellipsoids

    cs = [collect(Float64, E.c) for E in ellipsoids]
    nx = length(cs[1])

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])

        1 <= d <= nx || error("invalid periodic dimension: $d")
        p > 0 || error("invalid period (<= 0): $p")

        for k in 2:length(cs)
            delta_raw = cs[k][d] - cs[k - 1][d]
            delta_min = delta_raw - round(delta_raw / p) * p
            cs[k][d] = cs[k - 1][d] + delta_min
        end
    end

    out = UT.Ellipsoid[]
    for (i, E) in enumerate(ellipsoids)
        push!(out, UT.Ellipsoid(Matrix(E.P), cs[i]))
    end

    return out
end

"""
Wrap ellipsoid centers only (matrix P unchanged).
"""
function wrap_ellipsoid_centers(ellipsoids, periodic_dims, periodic_periods, periodic_start)
    isempty(ellipsoids) && return ellipsoids

    out = UT.Ellipsoid[]
    for E in ellipsoids
        c = collect(Float64, E.c)
        wrap_vector_periodic!(c, periodic_dims, periodic_periods, periodic_start)
        push!(out, UT.Ellipsoid(Matrix(E.P), c))
    end

    return out
end

"""
Attach unwrap metadata while preserving existing metadata payload.
"""
function augment_metadata_for_unwrap(metadata, was_unwrapped::Bool, wrap_jumps::Int)
    if metadata isa NamedTuple
        return merge(
            metadata,
            (;
                lmi_unwrapped = was_unwrapped,
                lmi_wrap_jumps = wrap_jumps,
            ),
        )
    end

    return (;
        original_metadata = metadata,
        lmi_unwrapped = was_unwrapped,
        lmi_wrap_jumps = wrap_jumps,
    )
end

"""
Return a certification-ready candidate trajectory (unwrap only when needed).
"""
function candidate_for_periodic_certification(
    candidate::OP.CandidateTrajectory,
    periodic_dims,
    periodic_periods,
)
    xs = collect(ST.enum_elems(candidate.x_traj))
    wrap_jumps = count_periodic_wrap_jumps(xs, periodic_dims, periodic_periods)

    if wrap_jumps == 0
        return candidate, false, wrap_jumps
    end

    xs_unwrapped = unwrap_periodic_state_list(xs, periodic_dims, periodic_periods)
    candidate_unwrapped = OP.CandidateTrajectory(
        ST.Trajectory(xs_unwrapped),
        candidate.u_traj;
        Ts = candidate.Ts,
        source = candidate.source,
        metadata = augment_metadata_for_unwrap(candidate.metadata, true, wrap_jumps),
    )

    return candidate_unwrapped, true, wrap_jumps
end
