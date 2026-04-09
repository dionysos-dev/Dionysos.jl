# =====================================================
#  TIME DISCRETIZATION FOR CONTINUOUS SYSTEMS
# =====================================================

"""
    runge_kutta4(dynamics, x, u, tstep, num_substeps::Int)

Perform one step of 4th-order Runge-Kutta time integration.

# Arguments
- `dynamics`: Function `f(x, u)` for the continuous dynamics
- `x`: Current state
- `u`: Control input
- `tstep`: Total time step
- `num_substeps`: Number of substeps for integration

# Returns
- Updated state after time `tstep`
"""
function runge_kutta4(dynamics, x, u, tstep, num_substeps::Int)
    τ = tstep / num_substeps
    for _ in 1:num_substeps
        k1 = dynamics(x, u)
        k2 = dynamics(x + k1 * (τ / 2), u)
        k3 = dynamics(x + k2 * (τ / 2), u)
        k4 = dynamics(x + k3 * τ, u)
        x = x + (k1 + 2 * k2 + 2 * k3 + k4) * (τ / 6)
    end
    return x
end

function simulate_control_map(dynamics; num_substeps = 5)
    return (x, u, tstep) -> runge_kutta4(dynamics, x, u, tstep, num_substeps)
end

function discretize_control_map(dynamics, tstep::Float64; num_substeps = 5)
    return (x, u) -> runge_kutta4(dynamics, x, u, tstep, num_substeps)
end

"""
    discretize_continuous_system(system::MS.AbstractContinuousSystem, tstep::Float64; num_substeps = 5)

Convert a continuous-time control system to a discrete-time system.

# Arguments
- `system`: Continuous-time `AbstractContinuousSystem` from MathematicalSystems.jl
- `tstep`: Fixed time step for Euler discretization
- `num_substeps`: Number of RK4 substeps per time step

# Returns
- A `ConstrainedBlackBoxControlDiscreteSystem` with the same state/input constraints
"""
function discretize_continuous_system(
    system::MS.AbstractContinuousSystem,
    tstep::Float64;
    num_substeps = 5,
)
    discretized_dynamics =
        discretize_control_map(MS.mapping(system), tstep; num_substeps = num_substeps)

    return MS.ConstrainedBlackBoxControlDiscreteSystem(
        discretized_dynamics,
        MS.statedim(system),
        MS.inputdim(system),
        MS.stateset(system),
        MS.inputset(system),
    )
end
