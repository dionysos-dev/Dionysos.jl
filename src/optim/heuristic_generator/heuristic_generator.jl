export AbstractHeuristicGenerator,
       set_problem!,
       generate!,
       get_trajectory,
       get_success,
       get_solve_time

abstract type AbstractHeuristicGenerator end

# Attach / update the concrete problem the generator must solve.
function set_problem!(gen::AbstractHeuristicGenerator, concrete_problem)
    error("set_problem! not implemented for $(typeof(gen))")
end

# Run the heuristic and store results internally.
function generate!(gen::AbstractHeuristicGenerator)
    error("generate! not implemented for $(typeof(gen))")
end

# Return the candidate trajectory (your CandidateTrajectory struct or any agreed type).
function get_trajectory(gen::AbstractHeuristicGenerator)
    error("get_trajectory not implemented for $(typeof(gen))")
end

# Optional common getters
get_success(gen::AbstractHeuristicGenerator) = false
get_solve_time(gen::AbstractHeuristicGenerator) = NaN
