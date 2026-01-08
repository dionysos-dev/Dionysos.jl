export ContinuousTrajectory, ContinuousTrajectoryAttribute
export DiscreteTrajectory
export last_mode

"""
    DiscreteTrajectory{Q, TT}

`q_0` is the starting mode and `transitions` is a sequence of discrete transitions in the system.
"""
struct DiscreteTrajectory{Q, TT}
    q_0::Q
    transitions::Vector{TT}
end

function DiscreteTrajectory{TT}(q_0) where {TT}
    return DiscreteTrajectory(q_0, TT[])
end

function last_mode(system, traj::DiscreteTrajectory)
    if isempty(traj.transitions)
        return traj.q_0
    else
        return HybridSystems.target(system, traj.transitions[end])
    end
end
Base.length(traj::DiscreteTrajectory) = length(traj.transitions)
function append(traj::DiscreteTrajectory, t)
    return DiscreteTrajectory(traj.q_0, [traj.transitions; t])
end

"""
    ContinuousTrajectory{T, XVT<:AbstractVector{T}, UVT<:AbstractVector{T}}

`x` is a sequence of points in the state space and `u` is a sequence of points in the input space.
"""
struct ContinuousTrajectory{T, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}
    x::Vector{XVT}
    u::Vector{UVT}
end

struct ContinuousTrajectoryAttribute <: MOI.AbstractModelAttribute end

"""
    HybridTrajectory{T, TT, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}

`discrete` is the discrete trajectory of type `DiscreteTrajectory` and  `continuous` is a the continuous trajectory of type `ContinuousTrajectory`.
"""
struct HybridTrajectory{T, TT, XVT <: AbstractVector{T}, UVT <: AbstractVector{T}}
    discrete::DiscreteTrajectory{TT}
    continuous::ContinuousTrajectory{T, XVT, UVT}
end

"""
    Trajectory{T}

provides the sequence of some elements of a trajectory.
"""
struct Trajectory{T}
    seq::Vector{T}
end

Base.length(traj::Trajectory) = length(traj.seq)
get_elem(traj::Trajectory, n::Int) = traj.seq[n]

@recipe function f(traj::Trajectory; dims = [1, 2])
    @series begin
        dims := dims
        UT.DrawTrajectory(traj.seq)
    end
end

function get_cost_trajectory(x_traj::Trajectory, u_traj::Trajectory, c)
    @assert length(x_traj) == length(u_traj) + 1

    cs = Float64[]
    total_cost = 0.0

    for i in 1:length(u_traj)
        ci = c(x_traj.seq[i], u_traj.seq[i])
        total_cost += ci
        push!(cs, ci)
    end

    return (c = Trajectory(cs), total_cost)
end

"""
    wrap_coord(x::SVector{N, T}, periodic_dims::SVector{P, Int}, periods::SVector{P, T}; start = zeros(SVector{P, T}))

Wraps the vector `x` into a periodic domain along dimensions specified in `periodic_dims`,
with period lengths `periods` and optional offset `start`.

# Arguments
- `x`: The coordinate vector to wrap.
- `periodic_dims`: Indices of the periodic dimensions.
- `periods`: Period lengths for the periodic dimensions.
- `start` (optional): Starting values of the periodic domains (defaults to `0.0`).

# Returns
A wrapped `SVector` where each periodic dimension is mapped to the interval `[start[i], start[i] + periods[i])`.
"""
function wrap_coord(
    x::SVector{N, T},
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T};
    start::SVector{P, T} = zeros(SVector{P, T}),
) where {N, P, T}
    return SVector{N, T}(ntuple(d -> begin
        i = findfirst(isequal(d), periodic_dims)
        if i === nothing
            x[d]
        else
            s = start[i]
            p = periods[i]
            mod(x[d] - s, p) + s
        end
    end, N))
end

function get_periodic_wrapper(
    periodic_dims::SVector{P, Int},
    periods::SVector{P, T};
    start::SVector{P, T} = zeros(SVector{P, T}),
) where {P, T}
    return x -> wrap_coord(x, periodic_dims, periods; start = start)
end

function get_closed_loop_trajectory(
    system,
    controller::MS.AbstractMap,
    x0,
    nstep::Integer;
    stopping = x -> false,
    wrap = identity,
    f_map_override = nothing,
)
    f = f_map_override == nothing ? MS.mapping(system) : f_map_override # expects (x,u)
    k = controller.h  # expects (x)

    x = wrap(x0)
    xs = Vector{typeof(x)}(undef, 0);
    push!(xs, x)
    us = Any[]  # could be typed if you know u type

    for _ in 1:nstep
        stopping(x) && break

        u = k(x)
        u === nothing && break

        x = wrap(f(x, u))

        push!(us, u)
        push!(xs, x)
    end

    return (x = Trajectory(xs), u = Trajectory(us))
end

# u = h(q, x)
# x^+ = f(x, u)
# q^+ = g(q, x) ou g(q, x^+) if update_on_next == true (it depends of the controller)
function get_closed_loop_trajectory(
    system,
    controller::MS.SystemWithOutput,
    x0,
    q0,
    nstep::Integer;
    meas = identity,
    stopping = x -> false,
    wrap = identity,
    update_on_next::Bool = false,
)
    f = MS.mapping(system)              # (x,u) -> xnext
    gc = MS.mapping(controller.s)        # (q,y) -> qnext
    hc = controller.outputmap.h # (q,y) -> u

    x = wrap(x0)
    q = q0

    xs = Vector{typeof(x)}(undef, 0);
    push!(xs, x)
    qs = Vector{typeof(q)}(undef, 0);
    push!(qs, q)
    us = Any[]

    for _ in 1:nstep
        stopping(x) && break

        y = meas(x)
        u = hc((q, y))
        u === nothing && break

        xnext = f(x, u)
        xnext = wrap(xnext)

        y_for_update = update_on_next ? meas(xnext) : y
        qnext = gc(q, y_for_update)

        x, q = xnext, qnext

        push!(us, u)
        push!(xs, x)
        push!(qs, q)
    end

    return (x = Trajectory(xs), u = Trajectory(us), q = Trajectory(qs))
end