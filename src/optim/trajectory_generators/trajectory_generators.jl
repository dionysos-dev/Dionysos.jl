export AbstractTrajectoryGenerator,
    set_problem!,
    set_seed_trajectory!,
    set_horizon!,
    set_stop_on_success!,
    generate!,
    get_trajectory,
    get_success,
    get_solve_time,
    get_seed,
    get_diagnostics

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

function get_success(gen::AbstractTrajectoryGenerator)
    return error("get_success not implemented for $(typeof(gen))")
end

function get_solve_time(gen::AbstractTrajectoryGenerator)
    return error("get_solve_time not implemented for $(typeof(gen))")
end

"""
    set_horizon!(gen::AbstractTrajectoryGenerator, nstep::Int)

Re-target the generator at an `nstep`-step horizon (the prefix re-planning loop
shortens the horizon to the failed prefix). Errors for generators without an
adjustable horizon.
"""
function set_horizon!(gen::AbstractTrajectoryGenerator, nstep::Int)
    return error("set_horizon! not implemented for $(typeof(gen))")
end

"""
    set_stop_on_success!(gen, flag::Bool) -> Union{Nothing, Bool}

Set whether the generator stops at its first success, returning the PREVIOUS
value so callers can restore it — or `nothing` when the generator has no such
switch (the default: not every generator optimizes past success).
"""
set_stop_on_success!(gen::AbstractTrajectoryGenerator, flag::Bool) = nothing

"""
    get_seed(gen::AbstractTrajectoryGenerator)

The generator's current seed trajectory, or `nothing`.
"""
function get_seed(gen::AbstractTrajectoryGenerator)
    return error("get_seed not implemented for $(typeof(gen))")
end

"""
    get_diagnostics(gen::AbstractTrajectoryGenerator)

Generator-specific diagnostics of the last `generate!` run, or `nothing`.
"""
function get_diagnostics(gen::AbstractTrajectoryGenerator)
    return error("get_diagnostics not implemented for $(typeof(gen))")
end

include("costs.jl")
include("rollout.jl")
include("optimizer_trajectory_generator.jl")
include("mppi/mppi.jl")
include("composite_trajectory_generator.jl")
