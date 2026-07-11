# --------------------------------------------------
#  CENTERED SIMULATION IMPLEMENTATION
# --------------------------------------------------

"""
    DiscreteTimeCenteredSimulation <: DiscreteTimeSystemUnderApproximation

A concrete underapproximation that simulates the evolution of the **center** of the input set under a discrete-time system.

This approximation is very conservative, returning a single propagated point from the center of the input set.

# Fields
- `system`: A constrained discrete-time control system (e.g., from [`MathematicalSystems.jl`](https://github.com/JuliaReach/MathematicalSystems.jl)).

# Underapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector) -> Vector{SVector}`
which returns a singleton list with the propagated center point.
"""
struct DiscreteTimeCenteredSimulation{S <: MS.ConstrainedBlackBoxControlDiscreteSystem} <:
       DiscreteTimeSystemUnderApproximation
    system::S
end

get_system(approx::DiscreteTimeCenteredSimulation) = approx.system

function get_under_approximation_map(approx::DiscreteTimeCenteredSimulation)
    system_map = get_system_map(approx)
    return (elem, u) -> begin
        x = LazySets.center(elem)
        Fx = system_map(x, u)
        return [Fx]  # Return a single point representing the propagated center
    end
end

"""
    ContinuousTimeCenteredSimulation <: ContinuousTimeSystemUnderApproximation

A concrete underapproximation of a continuous-time system using **center-point simulation**.

Simulates only the center of the state set under the system dynamics. Returns a single propagated point after integration over a time step.

# Fields
- `system`: A constrained continuous-time control system.

# Underapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector, tstep::Real) -> Vector{SVector}`
which returns a singleton list with the propagated center point.

# Notes
Use `discretize` to convert this approximation into a discrete-time approximation suitable for use in fixed-step abstraction pipelines.
"""
struct ContinuousTimeCenteredSimulation{
    S <: MS.ConstrainedBlackBoxControlContinuousSystem,
} <: ContinuousTimeSystemUnderApproximation
    system::S
end

get_system(approx::ContinuousTimeCenteredSimulation) = approx.system

function get_under_approximation_map(approx::ContinuousTimeCenteredSimulation)
    system_map = get_system_map(approx)
    return (elem, u, tstep) -> begin
        x = LazySets.center(elem)
        Fx = system_map(x, u, tstep)
        return [Fx]  # Return a single point representing the propagated center
    end
end

function discretize(
    approx::ContinuousTimeCenteredSimulation,
    tstep::Float64;
    num_substeps::Int = DEFAULT_NUM_SUBSTEPS,
)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = num_substeps)
    return DiscreteTimeCenteredSimulation(discretized_system)
end

# --------------------------------------------------
#  RANDOM SIMULATION IMPLEMENTATION
# --------------------------------------------------

"""
    DiscreteTimeRandomSimulation <: DiscreteTimeSystemUnderApproximation

A stochastic underapproximation of a discrete-time system using **random sampling**.

Propagates multiple randomly sampled points from the input set to provide a discrete underapproximation of reachable states.

# Fields
- `system`: The underlying discrete-time control system.
- `nsamples`: Number of samples to draw from the input set.

# Underapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector) -> Vector{SVector}`
which returns a list of propagated samples.
"""
struct DiscreteTimeRandomSimulation{S <: MS.ConstrainedBlackBoxControlDiscreteSystem} <:
       DiscreteTimeSystemUnderApproximation
    system::S
    nsamples::Int
end

get_system(approx::DiscreteTimeRandomSimulation) = approx.system

function get_under_approximation_map(approx::DiscreteTimeRandomSimulation)
    return (elem, u) -> begin
        samples = UT.samples(elem, approx.nsamples)
        return [get_system_map(approx)(x, u) for x in samples]
    end
end

"""
    ContinuousTimeRandomSimulation <: ContinuousTimeSystemUnderApproximation

A stochastic underapproximation for continuous-time systems using **random point sampling**.

Simulates multiple samples from the input set, over a fixed time step.

# Fields
- `system`: The underlying continuous-time control system.
- `nsamples`: Number of random samples.

# Underapproximation Map
Returns a function of the form:
    `f(rect::LazySets.AbstractHyperrectangle, u::SVector, tstep::Real) -> Vector{SVector}`
which returns a list of propagated samples.

# Notes
Use `discretize` to convert this approximation into a discrete-time approximation suitable for use in fixed-step abstraction pipelines.
"""
struct ContinuousTimeRandomSimulation{S <: MS.ConstrainedBlackBoxControlContinuousSystem} <:
       ContinuousTimeSystemUnderApproximation
    system::S
    nsamples::Int
end

get_system(approx::ContinuousTimeRandomSimulation) = approx.system

function get_under_approximation_map(approx::ContinuousTimeRandomSimulation)
    return (elem, u, tstep) -> begin
        samples = UT.samples(elem, approx.nsamples)
        return [get_system_map(approx)(x, u, tstep) for x in samples]
    end
end

function discretize(
    approx::ContinuousTimeRandomSimulation,
    tstep::Float64;
    num_substeps::Int = DEFAULT_NUM_SUBSTEPS,
)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = num_substeps)
    return DiscreteTimeRandomSimulation(discretized_system, approx.nsamples)
end
