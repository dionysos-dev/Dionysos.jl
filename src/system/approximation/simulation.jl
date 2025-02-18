# --------------------------------------------------
#  CENTERED SIMULATION IMPLEMENTATION
# --------------------------------------------------

struct DiscreteTimeCenteredSimulation <: DiscreteTimeSystemUnderApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlDiscreteSystem}
end

get_system(approx::DiscreteTimeCenteredSimulation) = approx.system

function get_under_approximation_map(approx::DiscreteTimeCenteredSimulation)
    system_map = get_system_map(approx)
    return (elem::Union{UT.HyperRectangle, UT.Ellipsoid}, u) -> begin
        x = UT.get_center(elem)
        Fx = system_map(x, u)
        return [Fx]  # Return a single point representing the propagated center
    end
end

struct ContinuousTimeCenteredSimulation <: ContinuousTimeSystemUnderApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlContinuousSystem}
end

get_system(approx::ContinuousTimeCenteredSimulation) = approx.system

function get_under_approximation_map(approx::ContinuousTimeCenteredSimulation)
    system_map = get_system_map(approx)
    return (elem::Union{UT.HyperRectangle, UT.Ellipsoid}, u, tstep) -> begin
        x = UT.get_center(elem)
        Fx = system_map(x, u, tstep)
        return [Fx]  # Return a single point representing the propagated center
    end
end

function discretize(approx::ContinuousTimeCenteredSimulation, tstep::Float64)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = 5)
    return DiscreteTimeCenteredSimulation(discretized_system)
end

# --------------------------------------------------
#  RANDOM SIMULATION IMPLEMENTATION
# --------------------------------------------------

struct DiscreteTimeRandomSimulation <: DiscreteTimeSystemUnderApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlDiscreteSystem}
    nsamples::Int
end

get_system(approx::DiscreteTimeRandomSimulation) = approx.system

function get_under_approximation_map(approx::DiscreteTimeRandomSimulation)
    return (elem::Union{UT.HyperRectangle, UT.Ellipsoid}, u) -> begin
        samples = UT.samples(elem, approx.nsamples)
        return [get_system_map(approx)(x, u) for x in samples]
    end
end

struct ContinuousTimeRandomSimulation <: ContinuousTimeSystemUnderApproximation
    system::Union{Nothing, MS.ConstrainedBlackBoxControlContinuousSystem}
    nsamples::Int
end

get_system(approx::ContinuousTimeRandomSimulation) = approx.system

function get_under_approximation_map(approx::ContinuousTimeRandomSimulation)
    return (elem::Union{UT.HyperRectangle, UT.Ellipsoid}, u, tstep) -> begin
        samples = UT.samples(elem, approx.nsamples)
        return [get_system_map(approx)(x, u, tstep) for x in samples]
    end
end

function discretize(approx::ContinuousTimeRandomSimulation, tstep::Float64)
    discretized_system =
        discretize_continuous_system(get_system(approx), tstep; num_substeps = 5)
    return DiscreteTimeRandomSimulation(discretized_system, approx.nsamples)
end
