# Closed-loop simulation directly on the abstract automaton.

function get_closed_loop_trajectory(
    autom::AbstractAutomatonList,
    controller::ST.AbstractDiscreteController,
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

    if trajectory_success(ST.Trajectory(qs))
        verbose &&
            @info "Closed-loop simulation stopped: trajectory success reached" step = 0 state =
                q
        return (x = ST.Trajectory(qs), u = ST.Trajectory(us))
    end

    for k in 1:nstep
        if stopping(q)
            verbose &&
                @info "Closed-loop simulation stopped: stopping condition reached" step = k state =
                    q
            break
        end

        u = ST.output_control(controller, nothing, q)

        if u === nothing
            verbose &&
                @warn "Closed-loop simulation stopped: controller returned nothing" step = k state =
                    q
            break
        end

        qnexts = post(autom, q, u)

        if qnexts === nothing || isempty(qnexts)
            verbose &&
                @warn "Closed-loop simulation stopped: no successor state" step = k state =
                    q input = u
            break
        end

        qnext = randomize_post ? rand(qnexts) : first(qnexts)

        push!(us, u)
        push!(qs, qnext)

        q = qnext

        if trajectory_success(ST.Trajectory(qs))
            verbose &&
                @info "Closed-loop simulation stopped: trajectory success reached" step = k state =
                    q
            break
        end
    end

    return (x = ST.Trajectory(qs), u = ST.Trajectory(us))
end
