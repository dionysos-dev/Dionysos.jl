export MPPIConfig, MPPIGenerator, get_seed

import Dionysos
import Random

const DI = Dionysos
const ST = DI.System

_identity_wrap_state(problem, x) = x

struct MPPIConfig{FX0, FDYN, FNOISE, FPROJ, FCOST, FSUCC, FWRAP, FPOST}
    Δt::Float64
    nstep::Int
    nsamples::Int
    niter::Int
    λ::Float64
    x0_provider::FX0
    discrete_dynamics::FDYN
    noise_sampler::FNOISE
    project_input::FPROJ
    trajectory_cost::FCOST
    success_fun::FSUCC
    wrap_state::FWRAP
    truncate_at_first_target::Bool
    postprocess_candidate::FPOST
end

function _default_mppi_postprocess(problem, cfg, cand)
    return cfg.truncate_at_first_target ?
           _truncate_candidate_at_first_target_hit(problem, cfg, cand) : cand
end

function MPPIConfig(
    Δt::Real,
    nstep::Integer,
    nsamples::Integer,
    niter::Integer,
    λ::Real,
    x0_provider,
    discrete_dynamics,
    noise_sampler,
    project_input,
    trajectory_cost,
    success_fun,
    wrap_state = _identity_wrap_state;
    truncate_at_first_target::Bool = true,
    postprocess_candidate = _default_mppi_postprocess,
)
    return MPPIConfig(
        Float64(Δt),
        Int(nstep),
        Int(nsamples),
        Int(niter),
        Float64(λ),
        x0_provider,
        discrete_dynamics,
        noise_sampler,
        project_input,
        trajectory_cost,
        success_fun,
        wrap_state,
        truncate_at_first_target,
        postprocess_candidate,
    )
end

mutable struct MPPIGenerator{
    P <: DI.Problem.ProblemType,
    C <: MPPIConfig,
    SG <: AbstractHeuristicGenerator,
} <: AbstractHeuristicGenerator
    problem::Union{Nothing, P}
    config::C
    seed_generator::SG
    seed::Union{Nothing, CandidateTrajectory}
    candidate::Union{Nothing, CandidateTrajectory}
    success::Bool
    solve_time_sec::Float64
    diagnostics::NamedTuple
end

function MPPIGenerator(
    problem::P,
    config::C,
    seed_generator::SG;
    diagnostics = (;),
) where {P <: DI.Problem.ProblemType, C <: MPPIConfig, SG <: AbstractHeuristicGenerator}
    return MPPIGenerator{P, C, SG}(
        problem,
        config,
        seed_generator,
        nothing,
        nothing,
        false,
        0.0,
        diagnostics,
    )
end

function MPPIGenerator(
    ::Nothing,
    config::C,
    seed_generator::SG;
    diagnostics = (;),
) where {C <: MPPIConfig, SG <: AbstractHeuristicGenerator}
    return MPPIGenerator{DI.Problem.ProblemType, C, SG}(
        nothing,
        config,
        seed_generator,
        nothing,
        nothing,
        false,
        0.0,
        diagnostics,
    )
end

function _pad_or_trim_controls(u_seq::AbstractVector, nstep::Int)
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
    problem::Union{Nothing, DI.Problem.ProblemType},
    cfg::MPPIConfig,
    x0,
    u_seq::AbstractVector,
)
    x = cfg.wrap_state(problem, x0)

    xs = [x]
    us = collect(u_seq)

    for (k, u) in enumerate(us)
        x = cfg.discrete_dynamics(problem, x, u, k, cfg.Δt)
        x = cfg.wrap_state(problem, x)
        push!(xs, x)
    end

    return CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = cfg.Δt,
        source = :mppi,
        metadata = (; nstep = cfg.nstep, nsamples = cfg.nsamples, niter = cfg.niter),
    )
end

function _truncate_candidate_at_first_target_hit(
    problem::DI.Problem.ProblemType,
    cfg::MPPIConfig,
    cand::Union{Nothing, CandidateTrajectory},
)
    cand === nothing && return nothing
    hasproperty(problem, :target_set) || return cand

    xs = collect(ST.enum_elems(cand.x_traj))
    hit_idx = findfirst(x -> (cfg.wrap_state(problem, x) ∈ problem.target_set), xs)

    if hit_idx === nothing || hit_idx <= 1 || hit_idx >= length(xs)
        return cand
    end

    us = collect(ST.enum_elems(cand.u_traj))

    return CandidateTrajectory(
        ST.Trajectory(xs[1:hit_idx]),
        ST.Trajectory(us[1:(hit_idx - 1)]);
        Ts = cand.Ts,
        source = cand.source,
        metadata = cand.metadata,
    )
end

function _sample_perturbed_controls(
    cfg::MPPIConfig,
    rng::Random.AbstractRNG,
    u_nom::AbstractVector,
)
    horizon = length(u_nom)
    eps_seq = [cfg.noise_sampler(rng, u_nom[k], k) for k in 1:horizon]
    u_roll = [cfg.project_input(u_nom[k] + eps_seq[k]) for k in 1:horizon]

    return eps_seq, u_roll
end

function _weighted_noise(eps_samples, weights, k::Int)
    δu = 0.0 * eps_samples[1][k]
    for s in eachindex(weights)
        δu += weights[s] * eps_samples[s][k]
    end
    return δu
end

function _mppi_update(
    problem::Union{Nothing, DI.Problem.ProblemType},
    cfg::MPPIConfig,
    x0,
    u_nom::AbstractVector,
    rng::Random.AbstractRNG,
)
    nsamples = cfg.nsamples
    horizon = length(u_nom)

    costs = Vector{Float64}(undef, nsamples)
    eps_seq, u_roll = _sample_perturbed_controls(cfg, rng, u_nom)
    cand = _rollout_controls(problem, cfg, x0, u_roll)
    costs[1] = Float64(cfg.trajectory_cost(problem, cand))
    eps_samples = [eps_seq]

    for s in 2:nsamples
        eps_seq, u_roll = _sample_perturbed_controls(cfg, rng, u_nom)
        cand = _rollout_controls(problem, cfg, x0, u_roll)
        costs[s] = Float64(cfg.trajectory_cost(problem, cand))
        push!(eps_samples, eps_seq)
    end

    β = minimum(costs)
    weights = exp.(-(costs .- β) ./ cfg.λ)
    weights ./= sum(weights)

    u_new = map(1:horizon) do k
        δu = _weighted_noise(eps_samples, weights, k)
        return cfg.project_input(u_nom[k] + δu)
    end

    return u_new, (; sample_best_cost = β)
end

function set_problem!(gen::MPPIGenerator, problem::DI.Problem.ProblemType)
    gen.problem = problem
    set_problem!(gen.seed_generator, problem)
    gen.seed = nothing
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    gen.diagnostics = (;)
    return gen
end

function generate!(gen::MPPIGenerator)
    cfg = gen.config

    @assert gen.problem !== nothing "Call set_problem!(gen, problem) first."
    @assert cfg.Δt > 0.0
    @assert cfg.nstep >= 1
    @assert cfg.nsamples >= 1
    @assert cfg.niter >= 1
    @assert cfg.λ > 0.0

    gen.seed = nothing
    gen.candidate = nothing
    gen.success = false
    gen.solve_time_sec = 0.0
    gen.diagnostics = (;)

    t0 = time()

    generate!(gen.seed_generator)
    seed = get_trajectory(gen.seed_generator)
    gen.seed = seed

    if seed === nothing
        gen.solve_time_sec = time() - t0
        gen.diagnostics = (; seed_available = false)
        return gen
    end

    problem = gen.problem
    x0 = cfg.x0_provider(problem)
    u_nom = collect(ST.enum_elems(seed.u_traj))
    u_nom = _pad_or_trim_controls(u_nom, cfg.nstep)
    u_nom = [cfg.project_input(u) for u in u_nom]

    best_cand = _rollout_controls(problem, cfg, x0, u_nom)
    seed_cost = Float64(cfg.trajectory_cost(problem, best_cand))
    best_cost = seed_cost

    rng = Random.default_rng()
    iterations_done = 0

    for it in 1:cfg.niter
        iterations_done = it

        u_new, _ = _mppi_update(problem, cfg, x0, u_nom, rng)

        cand = _rollout_controls(problem, cfg, x0, u_new)
        cost = Float64(cfg.trajectory_cost(problem, cand))

        u_nom = u_new
        if cost < best_cost
            best_cost = cost
            best_cand = cand
        end

        if cfg.success_fun(problem, best_cand)
            break
        end
    end

    final_cand = cfg.postprocess_candidate(problem, cfg, best_cand)

    gen.candidate = final_cand
    gen.success = cfg.success_fun(problem, final_cand)
    gen.solve_time_sec = time() - t0
    gen.diagnostics = (;
        seed_cost = seed_cost,
        best_cost = best_cost,
        final_cost = best_cost,
        iterations = iterations_done,
    )
    return gen
end

get_trajectory(gen::MPPIGenerator) = gen.candidate
get_seed(gen::MPPIGenerator) = gen.seed
get_success(gen::MPPIGenerator) = gen.success
get_solve_time(gen::MPPIGenerator) = gen.solve_time_sec
