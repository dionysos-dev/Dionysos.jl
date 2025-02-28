# --------------------------------------------------
#  GROWTH OVERAPPROXIMATION IMPLEMENTATION
# --------------------------------------------------

struct DiscreteTimeGrowthBound <: DiscreteTimeSystemOverApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlDiscreteSystem}
    growthbound_map::Function
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

struct ContinuousTimeGrowthBound <: ContinuousTimeSystemOverApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlContinuousSystem}
    growthbound_map::Function
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
function discretize(approx::ContinuousTimeGrowthBound, tstep::Float64)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = 5)
    discretized_growthbound_map = (r, u) -> approx.growthbound_map(r, u, tstep)
    return DiscreteTimeGrowthBound(discretized_system, discretized_growthbound_map)
end

function ContinuousTimeGrowthBound_from_jacobian_bound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    jacobian_bound;
    ngrowthbound = 5,
)
    modified_jacobian_bound = (r, u) -> jacobian_bound(u) * r
    growthbound_map =
        (r, u, tstep) -> RungeKutta4(modified_jacobian_bound, r, u, tstep, ngrowthbound)
    return ContinuousTimeGrowthBound(system, growthbound_map)
end

function ContinuousTimeGrowthBound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem;
    ngrowthbound = 5,
)
    jacobian_bound = compute_jacobian_bound(system)
    return ContinuousTimeGrowthBound_from_jacobian_bound(
        system,
        jacobian_bound;
        ngrowthbound = ngrowthbound,
    )
end

function compute_jacobian_bound(system::MS.ConstrainedBlackBoxControlContinuousSystem)
    # TODO: Implement the logic for computing the Jacobian bound.
end
