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

struct ClosedLoopTrajectory{XT, UT}
    x::Trajectory{XT}
    u::Trajectory{UT}
end

Base.length(traj::Trajectory) = length(traj.seq)
get_elem(traj::Trajectory, n::Int) = traj.seq[n]
enum_elems(traj::Trajectory) = traj.seq

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


function get_closed_loop_trajectory(
    autom::AbstractAutomatonList,
    controller::AbstractDiscreteController,
    q0::Int,
    nstep::Integer;
    stopping = q -> false,
    trajectory_success = traj -> false,
    randomize_post::Bool = false,
    verbose::Bool = false,
)
    q = q0

    qs = Int[]
    us = Int[]

    push!(qs, q)

    if trajectory_success(Trajectory(qs))
        verbose &&
            @info "Closed-loop simulation stopped: trajectory success reached" step = 0 state = q
        return (x = Trajectory(qs), u = Trajectory(us))
    end

    for k in 1:nstep
        if stopping(q)
            verbose &&
                @info "Closed-loop simulation stopped: stopping condition reached" step = k state = q
            break
        end

        u = output_control(controller, nothing, q)

        if u === nothing
            verbose &&
                @warn "Closed-loop simulation stopped: controller returned nothing" step = k state = q
            break
        end

        qnexts = post(autom, q, u)

        if qnexts === nothing || isempty(qnexts)
            verbose &&
                @warn "Closed-loop simulation stopped: no successor state" step = k state = q input = u
            break
        end

        qnext = randomize_post ? rand(qnexts) : first(qnexts)

        push!(us, u)
        push!(qs, qnext)

        q = qnext

        if trajectory_success(Trajectory(qs))
            verbose &&
                @info "Closed-loop simulation stopped: trajectory success reached" step = k state = q
            break
        end
    end

    return (x = Trajectory(qs), u = Trajectory(us))
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
    controller::AbstractController,
    x0,
    nstep::Integer;
    meas = identity,
    stopping = x -> false,
    trajectory_success = traj -> false,
    wrap = identity,
    update_on_next::Bool = false,
    f_map_override = nothing,
    verbose::Bool = false,
)
    f = f_map_override === nothing ? MS.mapping(system) : f_map_override

    x = wrap(x0)
    q = initial_state(controller)

    xs = Vector{typeof(x)}()
    push!(xs, x)

    us = Any[]

    if trajectory_success(Trajectory(xs))
        verbose &&
            @info "Closed-loop simulation stopped: trajectory success reached" step = 0 state = x
        return (x = Trajectory(xs), u = Trajectory(us))
    end

    if q === nothing
        for k in 1:nstep
            if stopping(x)
                verbose &&
                    @info "Closed-loop simulation stopped: stopping condition reached" step =
                        k state = x
                break
            end

            y = meas(x)
            u = output_control(controller, q, y)

            if u === nothing
                verbose &&
                    @warn "Closed-loop simulation stopped: controller returned nothing" step =
                        k state = x measurement = y controller_state = q
                break
            end

            xnext = wrap(f(x, u))

            if any(!isfinite, xnext)
                verbose &&
                    @warn "Closed-loop simulation stopped: non-finite next state" step = k state =
                        x input = u next_state = xnext
                break
            end

            x = xnext

            push!(us, u)
            push!(xs, x)

            if trajectory_success(Trajectory(xs))
                verbose &&
                    @info "Closed-loop simulation stopped: trajectory success reached" step =
                        k state = x
                break
            end
        end

        return (x = Trajectory(xs), u = Trajectory(us))
    else
        qs = Vector{typeof(q)}()
        push!(qs, q)

        for k in 1:nstep
            if stopping(x)
                verbose &&
                    @info "Closed-loop simulation stopped: stopping condition reached" step =
                        k state = x controller_state = q
                break
            end

            y = meas(x)
            u = output_control(controller, q, y)

            if u === nothing
                verbose &&
                    @warn "Closed-loop simulation stopped: controller returned nothing" step =
                        k state = x measurement = y controller_state = q
                break
            end

            xnext = wrap(f(x, u))

            if any(!isfinite, xnext)
                verbose &&
                    @warn "Closed-loop simulation stopped: non-finite next state" step = k state =
                        x input = u next_state = xnext controller_state = q
                break
            end

            y_for_update = update_on_next ? meas(xnext) : y
            qnext = update_state(controller, q, y_for_update)

            if qnext === nothing
                verbose &&
                    @warn "Closed-loop simulation stopped: controller state update returned nothing" step =
                        k state = x next_state = xnext measurement_for_update = y_for_update controller_state =
                        q
                break
            end

            x, q = xnext, qnext

            push!(us, u)
            push!(xs, x)
            push!(qs, q)

            if trajectory_success(Trajectory(xs))
                verbose &&
                    @info "Closed-loop simulation stopped: trajectory success reached" step =
                        k state = x controller_state = q
                break
            end
        end

        return (x = Trajectory(xs), u = Trajectory(us), q = Trajectory(qs))
    end
end
