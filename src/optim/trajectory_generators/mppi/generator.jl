# The MPPI trajectory generator: iterated open-loop MPPI with elitism, used as a
# one-shot trajectory optimizer on a discrete-time problem. Per iteration:
# pre-draw all noise sequentially (deterministic given the rng, independent of the
# thread schedule), roll out and score the samples in parallel, weight or refit,
# update the nominal, keep the best trajectory seen.

_identity_wrap_state(problem, x) = x

mutable struct TrajectoryGenerator{RNG, FNOISE, FPROJ, FCOST, FWRAP} <:
               AbstractTrajectoryGenerator
    # Inputs
    problem::Union{Nothing, PR.ProblemType} # system must be a MS.AbstractDiscreteSystem
    seed_generator::Union{Nothing, AbstractTrajectoryGenerator}
    seed_trajectory::Union{Nothing, ST.Trajectory}

    # Parameters
    rng::RNG
    nstep::Int
    nsamples::Int
    niter::Int
    λ::Float64
    noise::FNOISE          # GaussianMPPINoise or legacy (rng, u, k) -> ε closure
    project_input::FPROJ
    cost::FCOST            # CompositeCost / AbstractCostTerm or legacy (problem, traj) -> J
    wrap_state::FWRAP      # for periodic systems
    hard_constraint::Bool  # freeze rollouts leaving the state set (never truncate)
    violation_penalty::Float64
    ess_target::Union{Nothing, Float64}  # target ESS fraction; nothing = fixed λ
    α::Float64             # IS exploration parameter; α = 1 disables the cross term
    update_rule::Symbol    # :mppi or :cem
    elite_frac::Float64
    anneal::Float64        # per-iteration σ multiplier (structured noise only)
    antithetic::Bool       # mirrored ± noise pairs

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
    λ::Real = 1.0,
    noise = nothing,
    noise_sampler = nothing,
    project_input,
    cost = nothing,
    trajectory_cost = nothing,
    wrap_state = _identity_wrap_state,
    hard_constraint::Bool = false,
    violation_penalty::Real = 1e6,
    ess_target = nothing,
    α::Real = 1.0,
    update_rule::Symbol = :mppi,
    elite_frac::Real = 0.1,
    anneal::Real = 1.0,
    antithetic::Bool = false,
)
    noise_model = noise === nothing ? noise_sampler : noise
    noise_model === nothing &&
        error("Provide `noise` (GaussianMPPINoise) or `noise_sampler`.")
    cost_model = cost === nothing ? trajectory_cost : cost
    cost_model === nothing && error("Provide `cost` (cost terms) or `trajectory_cost`.")
    update_rule in (:mppi, :cem) || error("update_rule must be :mppi or :cem.")

    return TrajectoryGenerator{
        typeof(rng),
        typeof(noise_model),
        typeof(project_input),
        typeof(cost_model),
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
        noise_model,
        project_input,
        cost_model,
        wrap_state,
        hard_constraint,
        Float64(violation_penalty),
        ess_target === nothing ? nothing : Float64(ess_target),
        Float64(α),
        update_rule,
        Float64(elite_frac),
        Float64(anneal),
        antithetic,
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

# ------------------------------------------------------------
# Cost evaluation (function barrier over the two cost styles)
# ------------------------------------------------------------

_as_composite(c::CompositeCost) = c
_as_composite(t::AbstractCostTerm) = CompositeCost((t,))

function _evaluate_controls(gen, f, X, x0, us, wrapf)
    if gen.cost isa Union{CompositeCost, AbstractCostTerm}
        c, _ = rollout_cost(
            f,
            x0,
            us,
            wrapf,
            _as_composite(gen.cost),
            X,
            gen.hard_constraint,
            gen.violation_penalty,
        )
        return c
    end
    traj, frozen = rollout_trajectory(f, x0, us, wrapf, X, gen.hard_constraint)
    return Float64(gen.cost(gen.problem, traj)) + gen.violation_penalty * frozen
end

# ------------------------------------------------------------
# One update iteration
# ------------------------------------------------------------

function _mppi_update(gen::TrajectoryGenerator, f, X, x0, u_nom, wrapf, rng, σ_scale)
    ns = gen.nsamples
    horizon = length(u_nom)

    # Sequential noise draw: determinism is independent of the thread schedule.
    first_seq = [_draw(gen.noise, rng, u_nom[k], k, σ_scale) for k in 1:horizon]
    eps = Vector{typeof(first_seq)}(undef, ns)
    eps[1] = first_seq
    for s in 2:ns
        if gen.antithetic && iseven(s)
            eps[s] = [-e for e in eps[s - 1]]
        else
            eps[s] = [_draw(gen.noise, rng, u_nom[k], k, σ_scale) for k in 1:horizon]
        end
    end

    costs = Vector{Float64}(undef, ns)
    eff = Vector{typeof(first_seq)}(undef, ns)
    Threads.@threads for s in 1:ns
        u_roll = [gen.project_input(u_nom[k] + eps[s][k]) for k in 1:horizon]
        eff[s] = [u_roll[k] - u_nom[k] for k in 1:horizon]
        costs[s] = _evaluate_controls(gen, f, X, x0, u_roll, wrapf)
    end

    # Importance-sampling correction (structured noise, α < 1): γ·Σₖ ūₖᵀΣ⁻¹εₖ.
    if gen.α < 1.0 && _supports_is_term(gen.noise)
        γ = gen.λ * (1.0 - gen.α)
        for s in 1:ns
            costs[s] += γ * sum(_is_dot(gen.noise, u_nom[k], eff[s][k]) for k in 1:horizon)
        end
    end

    λ_used = gen.ess_target === nothing ? gen.λ : _ess_lambda(costs, gen.ess_target * ns)
    weights, ess, entropy = _softmin_weights(costs, λ_used)

    u_new =
        gen.update_rule === :cem ?
        _cem_new_controls(u_nom, eff, costs, gen.elite_frac, gen.project_input) :
        _mppi_new_controls(u_nom, eff, weights, gen.project_input)

    return u_new,
    (; sample_best_cost = minimum(costs), ess = ess, entropy = entropy, λ_used = λ_used)
end

# ------------------------------------------------------------
# The elitist outer loop
# ------------------------------------------------------------

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

    f = MS.mapping(problem.system)
    X = MS.stateset(problem.system)
    wrapf = x -> gen.wrap_state(problem, x)
    x0 = first(ST.states(seed))

    u_nom = collect(ST.inputs(seed))
    u_nom = _pad_or_trim_controls(u_nom, gen.nstep)
    u_nom = [gen.project_input(u) for u in u_nom]

    best_traj, _ = rollout_trajectory(f, x0, u_nom, wrapf, X, gen.hard_constraint)
    seed_cost = _evaluate_controls(gen, f, X, x0, u_nom, wrapf)
    best_cost = seed_cost

    iterations_done = 0
    last_info = (; sample_best_cost = NaN, ess = NaN, entropy = NaN, λ_used = gen.λ)

    for it in 1:gen.niter
        iterations_done = it
        σ_scale = gen.anneal^(it - 1)

        u_new, info = _mppi_update(gen, f, X, x0, u_nom, wrapf, gen.rng, σ_scale)
        last_info = info

        cand_cost = _evaluate_controls(gen, f, X, x0, u_new, wrapf)

        u_nom = u_new

        if cand_cost < best_cost
            best_cost = cand_cost
            best_traj, _ = rollout_trajectory(f, x0, u_new, wrapf, X, gen.hard_constraint)
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
        last_sample_best_cost = last_info.sample_best_cost,
        ess = last_info.ess,
        weight_entropy = last_info.entropy,
        λ_used = last_info.λ_used,
        iterations = iterations_done,
    )

    return gen
end

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

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
    @assert 0.0 < gen.α <= 1.0
    @assert 0.0 < gen.elite_frac <= 1.0
    @assert gen.anneal > 0.0
end
