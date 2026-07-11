# --------------------------------------------------
#  REACHABLE-SET APPROXIMATIONS OF THE DYNAMICS (UNDER AND OVER)
# --------------------------------------------------

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

get_system(approx::SystemApproximation) =
    error("implement `get_system` for $(typeof(approx))")
is_over_approximation(approx::SystemApproximation) =
    error("implement `is_over_approximation` for $(typeof(approx))")

is_continuous_time(approx::SystemApproximation) =
    isa(approx, ContinuousTimeSystemApproximation)

get_system_map(approx::DiscreteTimeSystemApproximation) = MS.mapping(get_system(approx))
get_system_map(
    approx::ContinuousTimeSystemApproximation;
    num_substeps::Int = DEFAULT_NUM_SUBSTEPS,
) = simulate_control_map(MS.mapping(get_system(approx)); num_substeps = num_substeps)

discretize(approx::ContinuousTimeSystemApproximation, tstep::Float64; kwargs...) =
    error("implement `discretize` for $(typeof(approx))")

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
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector{M,T}) -> SVector{N,T}[]`
"""
get_under_approximation_map(approx::DiscreteTimeSystemUnderApproximation) =
    error("implement `get_under_approximation_map` for $(typeof(approx))")

"""
    get_under_approximation_map(approx::ContinuousTimeSystemUnderApproximation) -> Function

Returns a function that computes the underapproximation (list of points) of the system's evolution:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector{M,T}, tstep::T) -> SVector{N,T}[]`
"""
get_under_approximation_map(approx::ContinuousTimeSystemUnderApproximation) =
    error("implement `get_under_approximation_map` for $(typeof(approx))")

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
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector{M,T}) -> LazySets.Hyperrectangle`
"""
get_over_approximation_map(approx::DiscreteTimeSystemOverApproximation) =
    error("implement `get_over_approximation_map` for $(typeof(approx))")

"""
    get_over_approximation_map(overApprox::ContinuousTimeSystemOverApproximation) -> Function

Returns a function that computes the overapproximation of the system's evolution:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector{M,T}, tstep::T) -> LazySets.Hyperrectangle`
"""
get_over_approximation_map(approx::ContinuousTimeSystemOverApproximation) =
    error("implement `get_over_approximation_map` for $(typeof(approx))")

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
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector) -> LazySets.Hyperrectangle`  
    which returns an overapproximated successor set.
"""
struct DiscreteTimeOverApproximationMap{
    S <: MS.ConstrainedBlackBoxControlDiscreteSystem,
    F,
} <: DiscreteTimeSystemOverApproximation
    system::S
    over_approximation_map::F
end
get_system(approx::DiscreteTimeOverApproximationMap) = approx.system
get_over_approximation_map(approx::DiscreteTimeOverApproximationMap) =
    approx.over_approximation_map

"""
    ContinuousTimeOverApproximationMap <: ContinuousTimeSystemOverApproximation

Concrete implementation of a continuous-time **overapproximation** of a control system.

This type stores a constrained continuous-time system and an overapproximation function that simulates or bounds the system’s behavior over a given time step.

# Fields
- `system`: The underlying `ConstrainedBlackBoxControlContinuousSystem` from `MathematicalSystems.jl`.
- `over_approximation_map`: A function of the form  
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector, tstep::Real) -> LazySets.Hyperrectangle`  
    which returns an overapproximated reachable set over the given time interval.

# Notes
Use `discretize` to convert this approximation into a discrete-time overapproximation suitable for use in fixed-step abstraction pipelines.
"""
struct ContinuousTimeOverApproximationMap{
    S <: MS.ConstrainedBlackBoxControlContinuousSystem,
    F,
} <: ContinuousTimeSystemOverApproximation
    system::S
    over_approximation_map::F
end
get_system(approx::ContinuousTimeOverApproximationMap) = approx.system
get_over_approximation_map(approx::ContinuousTimeOverApproximationMap) =
    approx.over_approximation_map
function discretize(
    approx::ContinuousTimeOverApproximationMap,
    tstep::Float64;
    num_substeps::Int = DEFAULT_NUM_SUBSTEPS,
)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = num_substeps)
    discrete_overapprox = (rect, u) -> get_over_approximation_map(approx)(rect, u, tstep)
    return DiscreteTimeOverApproximationMap(discretized_system, discrete_overapprox)
end

include("growth_bound.jl")
include("linearized.jl")
