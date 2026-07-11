# --------------------------------------------------
#  GROWTH OVERAPPROXIMATION IMPLEMENTATION
# --------------------------------------------------

"""
    DiscreteTimeGrowthBound <: DiscreteTimeSystemOverApproximation

A discrete-time overapproximation based on **growth bounds**.

Given a system and a `growthbound_map`, this approximation inflates the center trajectory by a radius that depends on the current state set's size and the input.

# Fields
- `system`: A `ConstrainedBlackBoxControlDiscreteSystem` from `MathematicalSystems.jl`.
- `growthbound_map`: A function  
    `f(radius::SVector, u::SVector) -> SVector`  
    that computes how uncertainty in state evolves under the system.

# Overapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector) -> LazySets.Hyperrectangle`
This function simulates the image of the center and inflates it using the computed growth bound.
"""
struct DiscreteTimeGrowthBound{S <: MS.ConstrainedBlackBoxControlDiscreteSystem, F} <:
       DiscreteTimeSystemOverApproximation
    system::S
    growthbound_map::F
end

get_system(approx::DiscreteTimeGrowthBound) = approx.system

input_cache(approx::DiscreteTimeGrowthBound, r, u) = approx.growthbound_map(r, u)

function reach_set(approx::DiscreteTimeGrowthBound, elem, u, Fr)
    Fx = get_system_map(approx)(LazySets.center(elem), u)
    return LazySets.Hyperrectangle(Fx, Fr)
end

function get_over_approximation_map(approx::DiscreteTimeGrowthBound)
    return (elem, u) -> begin
        r = LazySets.radius_hyperrectangle(elem)
        return reach_set(approx, elem, u, input_cache(approx, r, u))
    end
end

"""
    ContinuousTimeGrowthBound <: ContinuousTimeSystemOverApproximation

A continuous-time overapproximation based on **growth bounds** for reachable set propagation.

It estimates how uncertainty evolves through time using a `growthbound_map` which depends on the radius, input, and time step.

# Fields
- `system`: A `ConstrainedBlackBoxControlContinuousSystem` from `MathematicalSystems.jl`.
- `growthbound_map`: A function  
    `f(radius::SVector, u::SVector, tstep::Real) -> SVector`  
    that estimates how uncertainty grows over a time step.

# Overapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector, tstep::Real) -> LazySets.Hyperrectangle`
This function simulates the image of the center and inflates it using the computed growth bound.
"""
struct ContinuousTimeGrowthBound{S <: MS.ConstrainedBlackBoxControlContinuousSystem, F} <:
       ContinuousTimeSystemOverApproximation
    system::S
    growthbound_map::F
end

get_system(approx::ContinuousTimeGrowthBound) = approx.system
function get_over_approximation_map(approx::ContinuousTimeGrowthBound)
    return (rect, u, tstep) -> begin
        x = LazySets.center(rect)
        r = LazySets.radius_hyperrectangle(rect)
        Fx = get_system_map(approx)(x, u, tstep)
        Fr = approx.growthbound_map(r, u, tstep)
        return LazySets.Hyperrectangle(Fx, Fr)
    end
end
function discretize(
    approx::ContinuousTimeGrowthBound,
    tstep::Float64;
    num_substeps::Int = DEFAULT_NUM_SUBSTEPS,
)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = num_substeps)
    discretized_growthbound_map = (r, u) -> approx.growthbound_map(r, u, tstep)
    return DiscreteTimeGrowthBound(discretized_system, discretized_growthbound_map)
end

"""
    ContinuousTimeGrowthBound(system; jacobian_bound = nothing, ngrowthbound = DEFAULT_NUM_SUBSTEPS)

Build a continuous-time growth-bound overapproximation from a Jacobian bound: the
radius dynamics `ṙ = jacobian_bound(u) * r` are integrated with `ngrowthbound` RK4
substeps per time step. When `jacobian_bound` is not provided, it is derived from
the system via `compute_jacobian_bound` (requires an extension providing it).
"""
function ContinuousTimeGrowthBound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem;
    jacobian_bound = nothing,
    ngrowthbound::Int = DEFAULT_NUM_SUBSTEPS,
)
    if jacobian_bound === nothing
        jacobian_bound = compute_jacobian_bound(system)
    end
    modified_jacobian_bound = (r, u) -> jacobian_bound(u) * r
    growthbound_map =
        (r, u, tstep) -> runge_kutta4(modified_jacobian_bound, r, u, tstep, ngrowthbound)
    return ContinuousTimeGrowthBound(system, growthbound_map)
end

# Deprecated: use `ContinuousTimeGrowthBound(system; jacobian_bound, ngrowthbound)`.
function ContinuousTimeGrowthBound_from_jacobian_bound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    jacobian_bound;
    ngrowthbound::Int = DEFAULT_NUM_SUBSTEPS,
)
    return ContinuousTimeGrowthBound(
        system;
        jacobian_bound = jacobian_bound,
        ngrowthbound = ngrowthbound,
    )
end

# Untyped fallback so the Symbolics extension can add the typed method without
# overwriting it (its implementation traces the dynamics symbolically and bounds
# the Jacobian over X with interval arithmetic).
function compute_jacobian_bound(system)
    return error(
        "Automatic Jacobian-bound computation requires Symbolics.jl " *
        "(load it with `using Symbolics`), or provide the bound explicitly via " *
        "`ContinuousTimeGrowthBound(system; jacobian_bound = jacobian_bound)`.",
    )
end
