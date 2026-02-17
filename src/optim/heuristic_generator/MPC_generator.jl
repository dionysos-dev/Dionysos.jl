using Dionysos
const ST = Dionysos.System
using JuMP
using MathematicalSystems 
const MS = MathematicalSystems
using StaticArrays

mutable struct MPCGenerator <: AbstractHeuristicGenerator
    problem
    horizon::Int
    dt::Union{Nothing, Float64}  # only used for continuous-time systems
    optimizer

    trajectory::Union{Nothing, ST.Trajectory}
    success::Bool
    solve_time::Float64
end

function MPCGenerator(problem; horizon::Int, dt::Union{Nothing, Float64} = nothing, optimizer)
    horizon > 0 || error("Horizon must be positive.")
    dt !== nothing && dt <= 0 && error("Timestep must be positive.")

    return MPCGenerator(problem, horizon, dt, optimizer, nothing, false, NaN)
end

function MPCGenerator(; horizon::Int, dt::Union{Nothing, Float64} = nothing, optimizer)
    horizon > 0 || error("Horizon must be positive.")
    dt !== nothing && dt <= 0 && error("Timestep must be positive.")

    return MPCGenerator(nothing, horizon, dt, optimizer, nothing, false, NaN)
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
    horizon > 0 || error("Horizon must be positive.")
    gen.horizon = horizon
    return nothing
end

function set_timestep!(
    gen::MPCGenerator;
    dt::Float64,
)
    dt > 0 || error("Timestep must be positive.")
    gen.dt = dt
    return nothing
end

function is_discrete(system)
    return system isa MathematicalSystems.AbstractDiscreteSystem
end 

# Run the heuristic and store results internally.
function generate!(gen::MPCGenerator)

    # Ensure the problem is set 
    gen.problem === nothing && error("Problem not set for MPCGenerator. Call set_problem! first.")

    # Ensure parameters are set and valid
    gen.horizon === nothing && error("Horizon not set for MPCGenerator. Call set_horizon! first.")
    is_discrete(gen.problem.system) || (gen.dt === nothing && error("Timestep not set for MPCGenerator. Call set_timestep! first."))

    # reset old results
    gen.trajectory = nothing
    gen.success = false
    gen.solve_time = NaN

    t0 = time()
    solve_MPC!(gen)
    gen.solve_time = time() - t0

    println("MPC generation completed in $(gen.solve_time) seconds. Success: $(gen.success)")

    return nothing
end

function get_trajectory(gen::MPCGenerator)
    return gen.trajectory
end

get_success(gen::MPCGenerator) = gen.success
get_solve_time(gen::MPCGenerator) = gen.solve_time

function solve_MPC!(gen::MPCGenerator)

    # Extract necessary data
    problem = gen.problem
    N = gen.horizon
    dt = gen.dt
    target_set = problem.target_set

    system = problem.system

    # If the system is continuous, discretize it using the provided timestep
    !is_discrete(system) && (system = ST.discretize_continuous_system(system, dt))

    dynamics = MS.mapping(system)

    # Get initial state by sampling in the initial set 
    x0  = Dionysos.Utils.sample(problem.initial_set)
    nx = system.statedim
    nu = system.inputdim

    println("Solving MPC with horizon = $N and dt = $dt, nx = $nx, nu = $nu...")

    model = Model(gen.optimizer)

    @variable(model, x[1:nx, 1:N+1])  # State variables for each time step
    @variable(model, u[1:nu, 1:N])    # Control variables for each time step

    # Initial state constraint
    for i in 1:nx
         @constraint(model, x[i, 1] == x0[i])
    end

    # Dynamics constraints
    for k in 1:N 
        for i in 1:nx
            @NLconstraint(model, x[i, k+1] == dynamics(x[:, k], u[:, k])[i])
        end
    end

    # Domain constraints
    for k in 1:N+1
        @constraint(model, -pi <= x[1,k] <= pi)
        @constraint(model, -5.0 <= x[2,k] <= 5.0)
    end

    # Input set constraints 
    for k in 1:N
        for i in 1:nu
             @constraint(model, -4.0 <= u[i, k] <= 4.0)
        end
    end

    # Target set constraint at final time: for now hard constraint
    @constraint(model, pi - 15.0 * pi / 180.0 <= x[1, N+1] <= pi + 15.0 * pi / 180.0) # have to check how we could automate this
    @constraint(model, -1.0 <= x[2, N+1] <= 1.0)

    # Objective : 0 function just checks feasibility, can be replaced with actual cost function if desired
    @objective(model, Min, 0)

    # Solve the optimization problem
    optimize!(model)

    # Check solution status and extract trajectory if successful
    if JuMP.termination_status(model) == MOI.OPTIMAL || JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
        gen.success = true
        println("MPC optimization successful. Extracting trajectory...")
        println("value of x at final time: ", [value(x[i, N+1]) for i in 1:nx])
        # Extract the trajectory
        gen.trajectory = ST.Trajectory{SVector{nx, Float64}}([
            SVector{nx}(JuMP.value.(x[:,k]))
            for k in 1:N+1
        ])
    else
        gen.success = false
        error("MPC optimization did not find a solution. Status: $(JuMP.termination_status(model))")
    end

    return nothing
end