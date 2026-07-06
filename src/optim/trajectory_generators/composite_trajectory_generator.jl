module CompositeTrajectoryGenerator

import Dionysos

import ..AbstractTrajectoryGenerator
import ..set_problem!
import ..generate!
import ..get_trajectory
import ..get_success
import ..get_solve_time

import ..MPPITrajectoryGenerator: set_seed_trajectory!

const DI = Dionysos
const PR = DI.Problem

export TrajectoryGenerator, get_seed, get_diagnostics

mutable struct TrajectoryGenerator{
    GSEED <: AbstractTrajectoryGenerator,
    GREFINE <: AbstractTrajectoryGenerator,
} <: AbstractTrajectoryGenerator
    # Inputs
    problem::Union{Nothing, PR.ProblemType}
    seed_generator::GSEED
    refinement_generator::GREFINE
    Δt::Float64

    # Parameters
    num_substeps::Int

    # Outputs
    discrete_problem::Union{Nothing, PR.ProblemType}
    seed_trajectory::Any
    trajectory::Any
    success::Bool
    solve_time::Float64
    diagnostics::NamedTuple
end

function TrajectoryGenerator(
    seed_generator::GSEED,
    refinement_generator::GREFINE;
    Δt::Real,
    num_substeps::Integer = 5,
) where {GSEED <: AbstractTrajectoryGenerator, GREFINE <: AbstractTrajectoryGenerator}
    return TrajectoryGenerator{GSEED, GREFINE}(
        nothing,
        seed_generator,
        refinement_generator,
        Float64(Δt),
        Int(num_substeps),
        nothing,
        nothing,
        nothing,
        false,
        NaN,
        (;),
    )
end

function set_problem!(gen::TrajectoryGenerator, problem::PR.ProblemType)
    gen.problem = problem

    gen.discrete_problem =
        PR.discretize_problem(problem, gen.Δt; num_substeps = gen.num_substeps)

    set_problem!(gen.seed_generator, problem)
    set_problem!(gen.refinement_generator, gen.discrete_problem)

    gen.seed_trajectory = nothing
    gen.trajectory = nothing
    gen.success = false
    gen.solve_time = NaN
    gen.diagnostics = (;)

    return gen
end

function generate!(gen::TrajectoryGenerator)
    gen.problem === nothing && error("No problem attached. Call set_problem! first.")
    gen.discrete_problem === nothing &&
        error("Discrete problem was not built. Call set_problem! first.")

    t0 = time()

    generate!(gen.seed_generator)

    seed = get_trajectory(gen.seed_generator)
    gen.seed_trajectory = seed

    if seed === nothing
        gen.solve_time = time() - t0
        gen.success = false
        gen.diagnostics = (;
            seed_available = false,
            seed_success = get_success(gen.seed_generator),
            refinement_success = false,
        )
        return gen
    end

    set_seed_trajectory!(gen.refinement_generator, seed)
    generate!(gen.refinement_generator)

    gen.trajectory = get_trajectory(gen.refinement_generator)
    gen.success = get_success(gen.refinement_generator)
    gen.solve_time = time() - t0

    gen.diagnostics = (;
        seed_available = true,
        seed_success = get_success(gen.seed_generator),
        refinement_success = get_success(gen.refinement_generator),
        seed_solve_time = get_solve_time(gen.seed_generator),
        refinement_solve_time = get_solve_time(gen.refinement_generator),
        total_solve_time = gen.solve_time,
    )

    return gen
end

get_trajectory(gen::TrajectoryGenerator) = gen.trajectory
get_seed(gen::TrajectoryGenerator) = gen.seed_trajectory
get_success(gen::TrajectoryGenerator) = gen.success
get_solve_time(gen::TrajectoryGenerator) = gen.solve_time
get_diagnostics(gen::TrajectoryGenerator) = gen.diagnostics

end
