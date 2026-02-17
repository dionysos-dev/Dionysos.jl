mutable struct MPCGenerator <: AbstractHeuristicGenerator
    problem
    horizon::Int
    dt::Float64

    trajectory::Union{Nothing, ST.Trajectory}
    success::Bool
    solve_time::Float64
end

function MPCGenerator(problem; horizon::Int, dt::Float64)
    return MPCGenerator(problem, horizon, dt, nothing, false, NaN)
end

function MPCGenerator(; horizon::Int, dt::Float64)
    return MPCGenerator(nothing, horizon, dt, nothing, false, NaN)
end

function set_problem!(
    gen::MPCGenerator, 
    concrete_problem,
)
    gen.problem = concrete_problem
    # reset possible old results
    gen.trajectory = nothing
    gen.success = false
    gen.solve_time = NaN
    return nothing
end

function set_horizon!(
    gen::MPCGenerator;
    horizon::Int,
)
    gen.horizon = horizon
    return nothing
end

function set_timestep!(
    gen::MPCGenerator;
    dt::Float64,
)
    gen.dt = dt
    return nothing
end

# Run the heuristic and store results internally.
function generate!(gen::MPCGenerator)

    # Ensure the problem is set 
    gen.problem === nothing && error("Problem not set for MPCGenerator. Call set_problem! first.")

    # Ensure parameters are set and valid
    gen.horizon === nothing && error("Horizon not set for MPCGenerator. Call set_horizon! first.")
    gen.dt === nothing && error("Timestep not set for MPCGenerator. Call set_timestep! first.")
    gen.horizon <= 0 && error("Horizon must be positive for MPCGenerator.")
    gen.dt <= 0 && error("Timestep must be positive for MPCGenerator.")

    # reset old results
    gen.trajectory = nothing
    gen.success = false
    gen.solve_time = NaN

    t0 = time()
    solve_MPC!(gen)
    gen.solve_time = time() - t0

    !gen.success && @warn("MPC generation failed for horizon $(gen.horizon) and timestep $(gen.dt).")

    return nothing
end

function get_trajectory(gen::MPCGenerator)
    return gen.trajectory
end

get_success(gen::MPCGenerator) = gen.success
get_solve_time(gen::MPCGenerator) = gen.solve_time

function solve_MPC!(gen::MPCGenerator)

    

    return error("MPC solving logic not implemented for $(typeof(gen))")
end