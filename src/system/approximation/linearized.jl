# --------------------------------------------------
#  LINEARIZED OVERAPPROXIMATION IMPLEMENTATION
# --------------------------------------------------

"""
    DiscreteTimeLinearized <: DiscreteTimeSystemOverApproximation

A discrete-time overapproximation based on **linearization** of the system dynamics.

This model approximates system behavior by propagating the linearized dynamics around the center of the state set and bounding the resulting error.

# Fields
- `system`: A `ConstrainedBlackBoxControlDiscreteSystem` from `MathematicalSystems.jl`.
- `linsys_map`: A function `(x, dx, u) -> (Fx, DFx)` returning the linearized next state `Fx` and its Jacobian `DFx` around perturbation `dx`.
- `error_map`: A function `(radius, u) -> Δ` returning a bound on the linearization error based on the set radius.

# Overapproximation Map
Returns a function of the form:
    `f(rect::HyperRectangle, u::SVector) -> HyperRectangle`
It evaluates the system at the center, adds linearized spread based on Jacobian, and inflates with the error bound.
"""
struct DiscreteTimeLinearized <: DiscreteTimeSystemOverApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlDiscreteSystem}
    linsys_map::Function
    error_map::Function
end

get_system(approx::DiscreteTimeLinearized) = approx.system

function get_over_approximation_map(approx::DiscreteTimeLinearized)
    return (rect::UT.HyperRectangle, u) -> begin
        x = UT.get_center(rect)
        r = UT.get_r(rect)
        e = norm(r, Inf)
        N = UT.get_dims(rect)

        _H_ = SMatrix{N, N}(I) .* r
        _ONE_ = ones(SVector{N})

        Fe = approx.error_map(e, u)
        Fr = r .+ Fe

        Fx, DFx = approx.linsys_map(x, _H_, u)

        rad = abs.(DFx) * _ONE_ .+ Fe
        return UT.HyperRectangle(Fx - rad, Fx + rad)
    end
end

"""
    ContinuousTimeLinearized <: ContinuousTimeSystemOverApproximation

A continuous-time overapproximation based on **Runge-Kutta linearization** of the system dynamics.

The method propagates both the nominal trajectory and its linearized sensitivity over a time step using a 4th-order Runge-Kutta scheme, while bounding the second-order remainder error.

# Fields
- `system`: A `ConstrainedBlackBoxControlContinuousSystem` from `MathematicalSystems.jl`.
- `linsys_map`: A function `(x, dx, u, tstep) -> (Fx, DFx)` simulating a linearized trajectory and its Jacobian.
- `error_map`: A function `(r, u, tstep) -> Δ` computing a bound on the nonlinearity-induced error over time.

# Overapproximation Map
Returns a function of the form:
    `f(rect::HyperRectangle, u::SVector, tstep::Real) -> HyperRectangle`
The result is a conservative reachable set from the center using linearization + second-order error correction.
"""
struct ContinuousTimeLinearized <: ContinuousTimeSystemOverApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlContinuousSystem}
    linsys_map::Function
    error_map::Function
end

get_system(approx::ContinuousTimeLinearized) = approx.system

function get_over_approximation_map(approx::ContinuousTimeLinearized)
    return (rect::UT.HyperRectangle, u, tstep) -> begin
        x = UT.get_center(rect)
        r = UT.get_r(rect)
        e = norm(r, Inf)
        N = UT.get_dims(rect)

        _H_ = SMatrix{N, N}(I) .* r
        _ONE_ = ones(SVector{N})

        Fe = approx.error_map(e, u, tstep)
        Fr = r .+ Fe

        Fx, DFx = approx.linsys_map(x, _H_, u, tstep)
        rad = abs.(DFx) * _ONE_ .+ Fe
        return UT.HyperRectangle(Fx - rad, Fx + rad)
    end
end

function discretize(approx::ContinuousTimeLinearized, tstep::Float64)
    discretized_system = discretize_continuous_system(get_system(approx), tstep)
    linsys_map = (x, dx, u) -> approx.linsys_map(x, dx, u, tstep)
    error_map = (r, u) -> approx.error_map(r, u, tstep)
    return DiscreteTimeLinearized(discretized_system, linsys_map, error_map)
end

function get_DiscreteTimeGrowthBound(approx::DiscreteTimeLinearized)
    return DiscreteTimeGrowthBound(
        get_system(approx),
        approx.linsys_map,
        error_map.error_map,
    )
end

function get_ContinuousTimeGrowthBound(approx::ContinuousTimeLinearized)
    return ContinuousTimeGrowthBound(
        get_system(approx),
        approx.linsys_map,
        error_map.error_map,
    )
end

## Constructors
function RungeKutta4Linearized(F, DF, x, dx, u, tstep, num_substeps::Int)
    τ = tstep / num_substeps
    for _ in 1:num_substeps
        k1, Dk1 = F(x, u), DF(x, u) * dx
        k2, Dk2 = F(x + k1 * (τ / 2), u), DF(x + k1 * (τ / 2), u) * (dx + Dk1 * (τ / 2))
        k3, Dk3 = F(x + k2 * (τ / 2), u), DF(x + k2 * (τ / 2), u) * (dx + Dk2 * (τ / 2))
        k4, Dk4 = F(x + k3 * τ, u), DF(x + k3 * τ, u) * (dx + Dk3 * τ)

        x += (k1 + 2 * k2 + 2 * k3 + k4) * (τ / 6)
        dx += (Dk1 + 2 * Dk2 + 2 * Dk3 + Dk4) * (τ / 6)
    end
    return (x, dx)
end

function BoundSecondOrder(a, b, tstep)
    return a ≈ 0.0 ? b * tstep : (b / a) * exp(a * tstep) * (exp(a * tstep) - 1.0) / 2.0
end

function ContinuousTimeLinearized(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    DF_sys,
    bound_DF,
    bound_DDF;
    num_substeps = 5,
)
    linsys_map =
        (x, dx, u, tstep) ->
            RungeKutta4Linearized(system.f, DF_sys, x, dx, u, tstep, num_substeps)
    error_map = (r, u, tstep) -> BoundSecondOrder(bound_DF(u), bound_DDF(u), tstep) * (r^2)
    return ContinuousTimeLinearized(system, linsys_map, error_map)
end
