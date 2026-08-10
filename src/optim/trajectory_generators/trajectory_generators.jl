export AbstractTrajectoryGenerator,
    set_problem!,
    set_seed_trajectory!,
    generate!,
    get_trajectory,
    get_success,
    get_solve_time

abstract type AbstractTrajectoryGenerator end

"""
    set_seed_trajectory!(gen::AbstractTrajectoryGenerator, traj)

Warm-start `gen` with a seed trajectory. Part of the generator interface so a
composite (seed → refine) chain works with *any* refiner, not one hardcoded type.
"""
function set_seed_trajectory!(gen::AbstractTrajectoryGenerator, traj)
    return error("set_seed_trajectory! not implemented for $(typeof(gen))")
end

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

include("costs.jl")
include("rollout.jl")
include("optimizer_trajectory_generator.jl")
include("mppi/mppi.jl")
include("composite_trajectory_generator.jl")
