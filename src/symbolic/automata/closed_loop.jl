# Closed-loop simulation directly on the abstract automaton: the shared
# simulation engine in `System`, with the automaton `post` map (plus successor
# choice) as the step function.

function get_closed_loop_trajectory(
    autom::AbstractAutomatonList,
    controller::ST.AbstractController,
    q0::Int,
    nstep::Integer;
    stopping = q -> false,
    trajectory_success = traj -> false,
    randomize_post::Bool = false,
    verbose::Bool = false,
)
    step = (q, u) -> begin
        qnexts = post(autom, q, u)
        (qnexts === nothing || isempty(qnexts)) && return nothing
        return randomize_post ? rand(qnexts) : first(qnexts)
    end

    return ST.get_closed_loop_trajectory(
        nothing,
        controller,
        q0,
        nstep;
        f_map_override = step,
        stopping = stopping,
        trajectory_success = trajectory_success,
        verbose = verbose,
    )
end
