# Closed-loop simulation of a system under a controller — one engine for both
# controller kinds (static controllers run with `mem = nothing` and skip the
# memory update) and, via `f_map_override`, for other evolution models (the
# automaton simulation in `Symbolic` routes through here).

"""
    get_closed_loop_trajectory(system, controller, x0, nstep; kwargs...) -> ClosedLoopTrajectory

Simulate the closed loop for at most `nstep` steps and return the visited
states/inputs (plus the controller-memory trajectory `q` for dynamic
controllers). The loop stops early when `stopping(x)` holds, when
`trajectory_success(x_traj)` holds, when the controller returns `nothing`, or
when the successor is `nothing`/non-finite.

Keyword arguments: `meas` (measurement map, default `identity`), `stopping`,
`trajectory_success`, `wrap` (state normalization, e.g. periodic wrapping),
`update_on_next` (feed the *next* measurement to `update_state`),
`f_map_override` (replace `MathematicalSystems.mapping(system)` as the step
function; may return `nothing` for "no successor"), `verbose`.
"""
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
    isdyn = controller_kind(controller) isa DynamicKind
    f = f_map_override === nothing ? MS.mapping(system) : f_map_override

    x = wrap(x0)
    mem = initial_state(controller)

    xs = Vector{typeof(x)}()
    push!(xs, x)

    us = nothing
    mems = isdyn ? Vector{typeof(mem)}([mem]) : nothing

    _result() = ClosedLoopTrajectory(
        Trajectory(xs),
        Trajectory(us === nothing ? Any[] : us),
        mems === nothing ? nothing : Trajectory(mems),
    )

    if trajectory_success(Trajectory(xs))
        verbose &&
            @info "Closed-loop simulation stopped: trajectory success reached" step = 0 state =
                x
        return _result()
    end

    for k in 1:nstep
        if stopping(x)
            verbose &&
                @info "Closed-loop simulation stopped: stopping condition reached" step = k state =
                    x controller_state = mem
            break
        end

        y = meas(x)
        u = output_control(controller, mem, y)

        if u === nothing
            verbose &&
                @warn "Closed-loop simulation stopped: controller returned nothing" step = k state =
                    x measurement = y controller_state = mem
            break
        end

        xnext_raw = f(x, u)

        if xnext_raw === nothing
            verbose &&
                @warn "Closed-loop simulation stopped: no successor state" step = k state =
                    x input = u
            break
        end

        xnext = wrap(xnext_raw)

        if any(!isfinite, xnext)
            verbose &&
                @warn "Closed-loop simulation stopped: non-finite next state" step = k state =
                    x input = u next_state = xnext
            break
        end

        if isdyn
            y_for_update = update_on_next ? meas(xnext) : y
            memnext = update_state(controller, mem, y_for_update)

            if memnext === nothing
                verbose &&
                    @warn "Closed-loop simulation stopped: controller state update returned nothing" step =
                        k state = x next_state = xnext measurement_for_update = y_for_update controller_state =
                        mem
                break
            end

            mem = memnext
        end

        x = xnext

        us === nothing ? (us = [u]) : push!(us, u)
        push!(xs, x)
        mems === nothing || push!(mems, mem)

        if trajectory_success(Trajectory(xs))
            verbose &&
                @info "Closed-loop simulation stopped: trajectory success reached" step = k state =
                    x controller_state = mem
            break
        end
    end

    return _result()
end
