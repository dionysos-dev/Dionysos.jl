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
    `f(rect::HyperRectangle, u::SVector) -> HyperRectangle`
This function simulates the image of the center and inflates it using the computed growth bound.
"""
struct DiscreteTimeGrowthBound{S <: MS.ConstrainedBlackBoxControlDiscreteSystem, F} <:
       DiscreteTimeSystemOverApproximation
    system::S
    growthbound_map::F
end

get_system(approx::DiscreteTimeGrowthBound) = approx.system
function get_over_approximation_map(approx::DiscreteTimeGrowthBound)
    return (rect::UT.HyperRectangle, u) -> begin
        x = UT.get_center(rect)
        r = UT.get_r(rect)
        Fx = get_system_map(approx)(x, u)
        Fr = approx.growthbound_map(r, u)
        return UT.HyperRectangle(Fx - Fr, Fx + Fr)
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
    `f(rect::HyperRectangle, u::SVector, tstep::Real) -> HyperRectangle`
This function simulates the image of the center and inflates it using the computed growth bound.
"""
struct ContinuousTimeGrowthBound{S <: MS.ConstrainedBlackBoxControlContinuousSystem, F} <:
       ContinuousTimeSystemOverApproximation
    system::S
    growthbound_map::F
end

get_system(approx::ContinuousTimeGrowthBound) = approx.system
function get_over_approximation_map(approx::ContinuousTimeGrowthBound)
    return (rect::UT.HyperRectangle, u, tstep) -> begin
        x = UT.get_center(rect)
        r = UT.get_r(rect)
        Fx = get_system_map(approx)(x, u, tstep)
        Fr = approx.growthbound_map(r, u, tstep)
        return UT.HyperRectangle(Fx - Fr, Fx + Fr)
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

function ContinuousTimeGrowthBound_from_jacobian_bound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    jacobian_bound;
    ngrowthbound::Int = DEFAULT_NUM_SUBSTEPS,
)
    modified_jacobian_bound = (r, u) -> jacobian_bound(u) * r
    growthbound_map =
        (r, u, tstep) -> runge_kutta4(modified_jacobian_bound, r, u, tstep, ngrowthbound)
    return ContinuousTimeGrowthBound(system, growthbound_map)
end

function ContinuousTimeGrowthBound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem;
    ngrowthbound::Int = DEFAULT_NUM_SUBSTEPS,
)
    jacobian_bound = compute_jacobian_bound(system)
    return ContinuousTimeGrowthBound_from_jacobian_bound(
        system,
        jacobian_bound;
        ngrowthbound = ngrowthbound,
    )
end

function compute_jacobian_bound(system::MS.ConstrainedBlackBoxControlContinuousSystem)
    return error(
        "Automatic Jacobian-bound computation is not implemented. " *
        "Build the growth bound explicitly with a user-provided bound via " *
        "`ContinuousTimeGrowthBound_from_jacobian_bound(system, jacobian_bound)`.",
    )
end
