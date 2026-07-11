# Closed-loop simulation of a system under a controller.

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
