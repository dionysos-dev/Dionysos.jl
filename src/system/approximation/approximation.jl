# ------------------------------
#  UTILITIES
# ------------------------------

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

# --------------------------------------------------
#  SYSTEM APPROXIMATIONS (GENERAL)
# --------------------------------------------------

"""
    SystemApproximation

Abstract supertype for all system approximation types.  
"""
abstract type SystemApproximation end

"""
    DiscreteTimeSystemApproximation <: SystemApproximation

Abstract supertype for approximations of discrete-time systems.  
"""
abstract type DiscreteTimeSystemApproximation <: SystemApproximation end

"""
    ContinuousTimeSystemApproximation <: SystemApproximation

Abstract supertype for approximations of continuous-time systems.
"""
abstract type ContinuousTimeSystemApproximation <: SystemApproximation end

function get_system(approx::SystemApproximation) end
function is_over_approximation(approx::SystemApproximation) end

is_continuous_time(approx::SystemApproximation) =
    isa(approx, ContinuousTimeSystemApproximation)

get_system_map(approx::DiscreteTimeSystemApproximation) = MS.mapping(get_system(approx))
get_system_map(approx::ContinuousTimeSystemApproximation; num_substeps = 5) =
    simulate_control_map(MS.mapping(get_system(approx)); num_substeps = num_substeps)

function discretize(approx::ContinuousTimeSystemApproximation, tstep::Float64) end

# --------------------------------------------------
#  SYSTEM UNDERAPPROXIMATIONS
# --------------------------------------------------

"""
    DiscreteTimeSystemUnderApproximation <: DiscreteTimeSystemApproximation

An abstract type representing an **underapproximation** of a discrete-time system.  
"""
abstract type DiscreteTimeSystemUnderApproximation <: DiscreteTimeSystemApproximation end

"""
    ContinuousTimeSystemUnderApproximation <: ContinuousTimeSystemApproximation

An abstract type representing an **underapproximation** of a continuous-time system.  
"""
abstract type ContinuousTimeSystemUnderApproximation <: ContinuousTimeSystemApproximation end

is_over_approximation(approx::DiscreteTimeSystemUnderApproximation) = false
is_over_approximation(approx::ContinuousTimeSystemUnderApproximation) = false

"""
    get_under_approximation_map(approx::DiscreteTimeSystemUnderApproximation) -> Function

Returns a function that computes the underapproximation (list of points) of the system's evolution:
    `f(rect::UT.HyperRectangle{N,T}, u::SVector{M,T}) -> SVector{N,T}[]`
"""
function get_under_approximation_map(approx::DiscreteTimeSystemUnderApproximation) end

"""
    get_under_approximation_map(approx::ContinuousTimeSystemUnderApproximation) -> Function

Returns a function that computes the underapproximation (list of points) of the system's evolution:
    `f(rect::UT.HyperRectangle{N,T}, u::SVector{M,T}, tstep::T) -> SVector{N,T}[]`
"""
function get_under_approximation_map(approx::ContinuousTimeSystemUnderApproximation) end

include("simulation.jl")

# --------------------------------------------------
#  SYSTEM OVERAPPROXIMATIONS
# --------------------------------------------------

"""
    DiscreteTimeSystemOverApproximation <: DiscreteTimeSystemApproximation

An abstract type representing an **overapproximation** of a discrete-time system.  
"""
abstract type DiscreteTimeSystemOverApproximation <: DiscreteTimeSystemApproximation end

"""
    ContinuousTimeSystemOverApproximation <: ContinuousTimeSystemApproximation

An abstract type representing an **overapproximation** of a continuous-time system.
"""
abstract type ContinuousTimeSystemOverApproximation <: ContinuousTimeSystemApproximation end

is_over_approximation(approx::DiscreteTimeSystemOverApproximation) = true
is_over_approximation(approx::ContinuousTimeSystemOverApproximation) = true

"""
    get_over_approximation_map(approx::DiscreteTimeSystemOverApproximation) -> Function

Returns a function that computes the overapproximation of the system's evolution:
    `f(rect::UT.HyperRectangle{N,T}, u::SVector{M,T}) -> UT.HyperRectangle{N,T}`
"""
function get_over_approximation_map(approx::DiscreteTimeSystemOverApproximation) end

"""
    get_over_approximation_map(overApprox::ContinuousTimeSystemOverApproximation) -> Function

Returns a function that computes the overapproximation of the system's evolution:
    `f(rect::UT.HyperRectangle{N,T}, u::SVector{M,T}, tstep::T) -> UT.HyperRectangle{N,T}`
"""
function get_over_approximation_map(approx::ContinuousTimeSystemOverApproximation) end

function get_DiscreteTimeOverApproximationMap(approx::DiscreteTimeSystemOverApproximation)
    return DiscreteTimeOverApproximationMap(
        get_system(approx),
        get_over_approximation_map(approx),
    )
end
function get_DiscreteTimeOverApproximationMap(
    approx::ContinuousTimeSystemOverApproximation,
    tstep::Float64,
)
    discretized = discretize(approx, tstep)
    return get_DiscreteTimeOverApproximationMap(discretized)
end

# --------------------------------------------------
#  OVERAPPROXIMATION MAP IMPLEMENTATION
# --------------------------------------------------

"""
    DiscreteTimeOverApproximationMap <: DiscreteTimeSystemOverApproximation

Concrete implementation of a discrete-time **overapproximation** of a dynamical system.

This type wraps a constrained discrete-time system along with an overapproximation map that, given a set of states and a control input, returns a conservative reachable set.

# Fields
- `system`: The underlying `ConstrainedBlackBoxControlDiscreteSystem` from `MathematicalSystems.jl`.
- `over_approximation_map`: A function of the form  
    `f(rect::HyperRectangle, u::SVector) -> HyperRectangle`  
    which returns an overapproximated successor set.
"""
struct DiscreteTimeOverApproximationMap <: DiscreteTimeSystemOverApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlDiscreteSystem}
    over_approximation_map::Function
end
get_system(approx::DiscreteTimeOverApproximationMap) = approx.system
get_over_approximation_map(approx::DiscreteTimeOverApproximationMap) =
    approx.over_approximation_map

"""
    ContinuousTimeSystemOverApproximationMap <: ContinuousTimeSystemOverApproximation

Concrete implementation of a continuous-time **overapproximation** of a control system.

This type stores a constrained continuous-time system and an overapproximation function that simulates or bounds the system’s behavior over a given time step.

# Fields
- `system`: The underlying `ConstrainedBlackBoxControlContinuousSystem` from `MathematicalSystems.jl`.
- `over_approximation_map`: A function of the form  
    `f(rect::HyperRectangle, u::SVector, tstep::Real) -> HyperRectangle`  
    which returns an overapproximated reachable set over the given time interval.

# Notes
Use `discretize` to convert this approximation into a discrete-time overapproximation suitable for use in fixed-step abstraction pipelines.
"""
struct ContinuousTimeSystemOverApproximationMap <: ContinuousTimeSystemOverApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlContinuousSystem}
    over_approximation_map::Function
end
get_system(approx::ContinuousTimeSystemOverApproximationMap) = approx.system
get_over_approximation_map(approx::ContinuousTimeSystemOverApproximationMap) =
    approx.over_approximation_map
function discretize(approx::ContinuousTimeSystemOverApproximationMap, tstep::Float64)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = 5)
    discrete_overapprox = (rect, u) -> get_over_approximation_map(approx)(rect, u, tstep)
    return DiscreteTimeOverApproximationMap(discretized_system, discrete_overapprox)
end

include("growth.jl")
include("linearized.jl")
