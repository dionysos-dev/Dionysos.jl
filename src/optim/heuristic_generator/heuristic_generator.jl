include("candidate_trajectory.jl")

export AbstractHeuristicGenerator,
    set_problem!, generate!, get_trajectory, get_success, get_solve_time
export CandidateTrajectory, horizon, n_states
export CenteredAbstractionConfig, CenteredAbstractionGenerator
export MPPIConfig, MPPIGenerator

abstract type AbstractHeuristicGenerator end

# Attach / update the concrete problem the generator must solve.
function set_problem!(gen::AbstractHeuristicGenerator, concrete_problem)
    return error("set_problem! not implemented for $(typeof(gen))")
end

# Run the heuristic and store results internally.
function generate!(gen::AbstractHeuristicGenerator)
    return error("generate! not implemented for $(typeof(gen))")
end

# Return the candidate trajectory (your CandidateTrajectory struct or any agreed type).
function get_trajectory(gen::AbstractHeuristicGenerator)
    return error("get_trajectory not implemented for $(typeof(gen))")
end

# Optional common getters
get_success(gen::AbstractHeuristicGenerator) = false
get_solve_time(gen::AbstractHeuristicGenerator) = NaN

include("centered_abstraction_generator.jl")
include("mppi_generator.jl")
