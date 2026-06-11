import Dionysos
import StaticArrays: SVector

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const MP = DI.Mapping
const OP = DI.Optim
const SC = OP.SymbolicCertifier

const DEFAULT_PERIODIC_DIMS = SVector(3, 4)
const DEFAULT_PERIODIC_PERIODS = SVector(2pi, 2pi)
const DEFAULT_PERIODIC_START = SVector(-pi, -pi)

function _check_periodic_data(periodic_dims, periodic_periods)
    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims and periodic_periods must have the same length")
    return nothing
end

function _check_periodic_data(periodic_dims, periodic_periods, periodic_start)
    _check_periodic_data(periodic_dims, periodic_periods)
    length(periodic_dims) == length(periodic_start) ||
        error("periodic_dims and periodic_start must have the same length")
    return nothing
end

function unwrap_periodic_state_list(state_list, periodic_dims, periodic_periods)
    isempty(state_list) && return state_list
    _check_periodic_data(periodic_dims, periodic_periods)

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

function wrap_vector_periodic!(
    x::AbstractVector{Float64},
    periodic_dims,
    periodic_periods,
    periodic_start,
)
    _check_periodic_data(periodic_dims, periodic_periods, periodic_start)

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        s = Float64(periodic_start[i])
        x[d] = mod(x[d] - s, p) + s
    end
    return x
end

function wrap_periodic_state_list(
    state_list,
    periodic_dims,
    periodic_periods,
    periodic_start,
)
    isempty(state_list) && return state_list

    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]

    for x in xs
        wrap_vector_periodic!(x, periodic_dims, periodic_periods, periodic_start)
    end

    return [SVector{nx, Float64}(x) for x in xs]
end

function count_periodic_wrap_jumps(state_list, periodic_dims, periodic_periods)
    length(state_list) <= 1 && return 0
    _check_periodic_data(periodic_dims, periodic_periods)

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

function unwrap_ellipsoid_centers(ellipsoids, periodic_dims, periodic_periods)
    isempty(ellipsoids) && return ellipsoids

    cs = [SVector{length(E.c), Float64}(E.c) for E in ellipsoids]
    cs_unwrapped = unwrap_periodic_state_list(cs, periodic_dims, periodic_periods)

    return [
        UT.Ellipsoid(Matrix(E.P), collect(cs_unwrapped[i])) for
        (i, E) in enumerate(ellipsoids)
    ]
end

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

function augment_metadata_for_unwrap(metadata, was_unwrapped::Bool, wrap_jumps::Int)
    payload = (; lmi_unwrapped = was_unwrapped, lmi_wrap_jumps = wrap_jumps)

    if metadata isa NamedTuple
        return merge(metadata, payload)
    end

    return (; original_metadata = metadata, payload...)
end

function unwrap_candidate_for_certification(
    candidate::OP.CandidateTrajectory;
    periodic_dims,
    periodic_periods,
)
    xs = collect(ST.enum_elems(candidate.x_traj))
    wrap_jumps = count_periodic_wrap_jumps(xs, periodic_dims, periodic_periods)

    wrap_jumps == 0 && return candidate

    xs_unwrapped = unwrap_periodic_state_list(xs, periodic_dims, periodic_periods)

    return OP.CandidateTrajectory(
        ST.Trajectory(xs_unwrapped),
        candidate.u_traj;
        Ts = candidate.Ts,
        source = candidate.source,
        metadata = augment_metadata_for_unwrap(candidate.metadata, true, wrap_jumps),
    )
end

function build_periodic_certification_preprocessor(; periodic_dims, periodic_periods)
    return candidate -> unwrap_candidate_for_certification(
        candidate;
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
    )
end
