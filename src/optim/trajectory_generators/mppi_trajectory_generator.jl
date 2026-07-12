module MPPITrajectoryGenerator

import MathematicalSystems as MS
import Dionysos

import ..AbstractTrajectoryGenerator
import ..set_problem!
import ..generate!
import ..get_trajectory
import ..get_success
import ..get_solve_time

const DI = Dionysos
const ST = DI.System
const PR = DI.Problem

_identity_wrap_state(problem, x) = x

mutable struct TrajectoryGenerator{RNG, FNOISE, FPROJ, FCOST, FWRAP} <:
               AbstractTrajectoryGenerator
    # Inputs
    problem::Union{Nothing, PR.ProblemType} # problems with system being a MS.AbstractDiscreteSystem
    seed_generator::Union{Nothing, AbstractTrajectoryGenerator}
    seed_trajectory::Union{Nothing, ST.Trajectory}

    # Parameters
    rng::RNG
    nstep::Int
    nsamples::Int
    niter::Int
    λ::Float64

    noise_sampler::FNOISE  # sample noise for each control input at each time step
    project_input::FPROJ   # ensure that every perturbed control remains admissible
    trajectory_cost::FCOST # cost function for evaluating trajectories
    wrap_state::FWRAP      # for periodic systems
    hard_constraint::Bool  # whether to enforce hard state constraints

    # Outputs
    trajectory::Union{Nothing, ST.Trajectory}
    success::Bool
    solve_time_sec::Float64
    diagnostics::NamedTuple
end

function TrajectoryGenerator(;
    rng,
    seed_generator = nothing,
    seed_trajectory = nothing,
    nstep::Integer,
    nsamples::Integer,
    niter::Integer,
    λ::Real,
    noise_sampler,
    project_input,
    trajectory_cost,
    wrap_state = _identity_wrap_state,
    hard_constraint::Bool = true,
)
    return TrajectoryGenerator{
        typeof(rng),
        typeof(noise_sampler),
        typeof(project_input),
        typeof(trajectory_cost),
        typeof(wrap_state),
    }(
        nothing,
        seed_generator,
        seed_trajectory,
        rng,
        Int(nstep),
        Int(nsamples),
        Int(niter),
        Float64(λ),
        noise_sampler,
        project_input,
        trajectory_cost,
        wrap_state,
        hard_constraint,
        nothing,
        false,
        NaN,
        (;),
    )
end

function set_problem!(gen::TrajectoryGenerator, problem::PR.ProblemType)
    gen.problem = problem

    if gen.seed_generator !== nothing
        set_problem!(gen.seed_generator, problem)
    end

    gen.trajectory = nothing
    gen.success = false
    gen.solve_time_sec = NaN
    gen.diagnostics = (;)

    return gen
end

function set_seed_trajectory!(gen::TrajectoryGenerator, traj::ST.Trajectory)
    gen.seed_trajectory = traj
    return gen
end

get_trajectory(gen::TrajectoryGenerator) = gen.trajectory
get_seed(gen::TrajectoryGenerator) = gen.seed_trajectory
get_success(gen::TrajectoryGenerator) = gen.success
get_solve_time(gen::TrajectoryGenerator) = gen.solve_time_sec
get_diagnostics(gen::TrajectoryGenerator) = gen.diagnostics

function generate!(gen::TrajectoryGenerator)
    _validate!(gen)

    problem = gen.problem

    gen.trajectory = nothing
    gen.success = false
    gen.solve_time_sec = NaN
    gen.diagnostics = (;)

    t0 = time()

    seed = _get_seed_trajectory(gen)
    gen.seed_trajectory = seed

    if seed === nothing
        gen.solve_time_sec = time() - t0
        gen.diagnostics = (; seed_available = false)
        return gen
    end

    isempty(ST.states(seed)) && error("Seed trajectory has no state sequence.")

    x0 = first(ST.states(seed))

    u_nom = collect(ST.inputs(seed))
    u_nom = _pad_or_trim_controls(u_nom, gen.nstep)
    u_nom = [gen.project_input(u) for u in u_nom]

    best_traj = _rollout_controls(gen, problem, x0, u_nom)
    seed_cost = Float64(gen.trajectory_cost(problem, best_traj))
    best_cost = seed_cost

    iterations_done = 0
    last_sample_best_cost = NaN

    for it in 1:gen.niter
        iterations_done = it

        u_new, info = _mppi_update(gen, problem, x0, u_nom, gen.rng)
        last_sample_best_cost = info.sample_best_cost

        cand = _rollout_controls(gen, problem, x0, u_new)
        cost = Float64(gen.trajectory_cost(problem, cand))

        u_nom = u_new

        if cost < best_cost
            best_cost = cost
            best_traj = cand
        end

        if PR.trajectory_success(problem, best_traj)
            break
        end
    end

    best_traj = truncate_at_target(problem, best_traj)
    gen.trajectory = best_traj
    gen.success = PR.trajectory_success(problem, best_traj)
    gen.solve_time_sec = time() - t0
    gen.diagnostics = (;
        seed_available = true,
        seed_cost = seed_cost,
        best_cost = best_cost,
        final_cost = best_cost,
        last_sample_best_cost = last_sample_best_cost,
        iterations = iterations_done,
    )

    return gen
end

function _pad_or_trim_controls(u_seq::AbstractVector, nstep::Int)
    isempty(u_seq) && error("Seed trajectory has no input sequence.")

    length(u_seq) == nstep && return collect(u_seq)
    length(u_seq) > nstep && return collect(u_seq[1:nstep])

    padded = collect(u_seq)
    last_u = padded[end]

    while length(padded) < nstep
        push!(padded, last_u)
    end

    return padded
end

function _rollout_controls(
    gen::TrajectoryGenerator,
    problem::PR.ProblemType,
    x0,
    u_seq::AbstractVector,
)
    discrete_time_system = problem.system
    f = MS.mapping(discrete_time_system)
    x = gen.wrap_state(problem, x0)

    xs = [x]
    us = collect(u_seq)

    for (_, u) in enumerate(us)
        x = f(x, u)
        x = gen.wrap_state(problem, x)
        if gen.hard_constraint && !(x ∈ MS.stateset(discrete_time_system))
            break
        end
        push!(xs, x)
    end

    return ST.Trajectory(xs; inputs = us)
end

function _sample_perturbed_controls(gen::TrajectoryGenerator, rng, u_nom::AbstractVector)
    horizon = length(u_nom)

    eps_seq = [gen.noise_sampler(rng, u_nom[k], k) for k in 1:horizon]

    u_roll = [gen.project_input(u_nom[k] + eps_seq[k]) for k in 1:horizon]

    return eps_seq, u_roll
end

function _weighted_noise(eps_samples, weights, k::Int)
    δu = zero(eps_samples[1][k])

    for s in eachindex(weights)
        δu = δu + weights[s] * eps_samples[s][k]
    end

    return δu
end

function _mppi_update(
    gen::TrajectoryGenerator,
    problem::PR.ProblemType,
    x0,
    u_nom::AbstractVector,
    rng,
)
    nsamples = gen.nsamples
    horizon = length(u_nom)

    costs = Vector{Float64}(undef, nsamples)
    eps_samples = Vector{Any}(undef, nsamples)

    for s in 1:nsamples
        eps_seq, u_roll = _sample_perturbed_controls(gen, rng, u_nom)
        cand = _rollout_controls(gen, problem, x0, u_roll)

        costs[s] = Float64(gen.trajectory_cost(problem, cand))
        eps_samples[s] = eps_seq
    end

    β = minimum(costs)

    weights = exp.(-(costs .- β) ./ gen.λ)

    Z = sum(weights)
    if !isfinite(Z) || Z <= 0
        weights .= 1 / length(weights)
    else
        weights ./= Z
    end

    u_new = map(1:horizon) do k
        δu = _weighted_noise(eps_samples, weights, k)
        return gen.project_input(u_nom[k] + δu)
    end

    return u_new, (; sample_best_cost = β)
end

function _get_seed_trajectory(gen::TrajectoryGenerator)
    if gen.seed_trajectory !== nothing
        return gen.seed_trajectory
    end

    if gen.seed_generator !== nothing
        generate!(gen.seed_generator)
        return get_trajectory(gen.seed_generator)
    end

    return nothing
end

function truncate_at_target(problem, traj::ST.Trajectory)
    xs = ST.states(traj)
    idx = findfirst(x -> x ∈ problem.target_set, xs)

    idx === nothing && return traj
    idx == length(xs) && return traj

    return ST.Trajectory(xs[1:idx]; inputs = ST.inputs(traj)[1:(idx - 1)])
end

function _validate!(gen::TrajectoryGenerator)
    @assert gen.problem !== nothing "Call set_problem!(gen, problem) first."
    @assert gen.problem.system isa MS.AbstractDiscreteSystem "Should be a discrete-time system."
    @assert gen.nstep >= 1
    @assert gen.nsamples >= 1
    @assert gen.niter >= 1
    @assert gen.λ > 0.0
end

end # module
