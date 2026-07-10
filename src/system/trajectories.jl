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

get_closed_loop_trajectory(
    system,
    controller::AbstractController,
    x0,
    nstep::Integer;
    kwargs...,
) = _closed_loop(controller_kind(controller), system, controller, x0, nstep; kwargs...)

function _closed_loop(
    ::StaticKind,
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

    xs = Vector{typeof(x)}()
    push!(xs, x)

    us = Any[]

    if trajectory_success(Trajectory(xs))
        verbose &&
            @info "Closed-loop simulation stopped: trajectory success reached" step = 0 state =
                x
        return (x = Trajectory(xs), u = Trajectory(us))
    end

    for k in 1:nstep
        if stopping(x)
            verbose &&
                @info "Closed-loop simulation stopped: stopping condition reached" step = k state =
                    x
            break
        end

        y = meas(x)
        u = output_control(controller, nothing, y)

        if u === nothing
            verbose &&
                @warn "Closed-loop simulation stopped: controller returned nothing" step = k state =
                    x measurement = y
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
                @info "Closed-loop simulation stopped: trajectory success reached" step = k state =
                    x
            break
        end
    end

    return (x = Trajectory(xs), u = Trajectory(us))
end

function _closed_loop(
    ::DynamicKind,
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
            @info "Closed-loop simulation stopped: trajectory success reached" step = 0 state =
                x
        return (x = Trajectory(xs), u = Trajectory(us))
    end

    qs = Vector{typeof(q)}()
    push!(qs, q)

    for k in 1:nstep
        if stopping(x)
            verbose &&
                @info "Closed-loop simulation stopped: stopping condition reached" step = k state =
                    x controller_state = q
            break
        end

        y = meas(x)
        u = output_control(controller, q, y)

        if u === nothing
            verbose &&
                @warn "Closed-loop simulation stopped: controller returned nothing" step = k state =
                    x measurement = y controller_state = q
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
                @info "Closed-loop simulation stopped: trajectory success reached" step = k state =
                    x controller_state = q
            break
        end
    end

    return (x = Trajectory(xs), u = Trajectory(us), q = Trajectory(qs))
end
