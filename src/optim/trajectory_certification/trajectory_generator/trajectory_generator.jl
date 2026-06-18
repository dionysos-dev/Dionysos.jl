export AbstractTrajectoryGenerator,
    set_problem!,
    generate!,
    get_trajectory,
    get_success,
    get_solve_time

abstract type AbstractTrajectoryGenerator end

function set_problem!(gen::AbstractTrajectoryGenerator, concrete_problem)
    return error("set_problem! not implemented for $(typeof(gen))")
end

function generate!(gen::AbstractTrajectoryGenerator)
    return error("generate! not implemented for $(typeof(gen))")
end

function get_trajectory(gen::AbstractTrajectoryGenerator)
    return error("get_trajectory not implemented for $(typeof(gen))")
end

get_success(gen::AbstractTrajectoryGenerator) = false
get_solve_time(gen::AbstractTrajectoryGenerator) = NaN

include("optimizer_trajectory_generator.jl")
include("mppi_trajectory_generator.jl")
include("composite_trajectory_generator.jl")