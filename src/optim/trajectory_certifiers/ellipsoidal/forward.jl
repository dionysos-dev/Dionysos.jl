# Forward certification (plan.md §4.4): the tube is *propagated* from a given entry
# ellipsoid, one `ST.solve_transition_forward` SDP per step. The source of every
# step is data, so the linearization box is known exactly before solving — the
# x-side consistency gate holds by construction and no adaptive boxes exist here.
# Failure is late but informative (the contraction profile shows where the tube
# inflates); the terminal gate `E_{K+1} ⊆ target_set` closes the specification.

"""
    ForwardOptions(; kwargs...)

Options of the forward certification chain.

- `target_mode` — `:free` (target shape `Q₂` is a decision variable, trace
  objective, conditioning sandwich `q_min`/`q_max`) or `:fixed` (shape follows the
  entry ellipsoid's, only the scale `α` is free — its per-step value is the
  contraction profile);
- `α_max` — fail-fast guard: a step whose tube scale exceeds this (relative to the
  entry shape) is rejected (`Inf` disables);
- `entry_shape` — LazySets shape matrix of the entry ellipsoid `E₁`, or `nothing`
  to circumscribe the problem's initial set (centered at the trajectory start);
- `maxδu`, `λ`, `transition_cost`, `linearization_δu` as in [`ChainOptions`](@ref);
  `linearization_δx_margin ≥ 1` inflates the (known) state box handed to the
  Hessian bound, buying u-side slack;
- gates: `r_min`, `check_state_domain`, `check_terminal` as in [`ChainOptions`](@ref).
"""
Base.@kwdef mutable struct ForwardOptions
    target_mode::Symbol = :free
    α_max::Float64 = Inf
    q_min::Float64 = 1e-9
    q_max::Float64 = 1e9
    entry_shape::Union{Nothing, Matrix{Float64}} = nothing
    maxδu::Float64 = 0.5
    λ::Float64 = 0.01
    transition_cost::Union{Nothing, UT.QuadraticStateControlFunction, Matrix{Float64}} =
        nothing
    linearization_δu::Vector{Float64} = Float64[]
    linearization_δx_margin::Float64 = 1.1
    r_min::Float64 = 0.0
    check_state_domain::Bool = true
    check_terminal::Bool = true
    remainder_model::Symbol = :vertices
end

# Circumscribed ellipsoid of the initial set, centered at the trajectory start:
# for a box with vertices vⱼ, the diagonal shape with semi-axes √n·maxⱼ|vⱼᵢ − cᵢ|
# contains every vertex (Cauchy–Schwarz), hence the box.
function _entry_ellipsoid(problem, x1, opts::ForwardOptions)
    if opts.entry_shape !== nothing
        return LazySets.Ellipsoid(collect(float.(x1)), Matrix{Float64}(opts.entry_shape)),
        nothing
    end
    verts = try
        LazySets.vertices_list(problem.initial_set)
    catch
        return nothing,
        "cannot circumscribe the initial set (no vertex list) — set entry_shape"
    end
    isempty(verts) && return nothing, "empty initial set"
    n = length(x1)
    h = [maximum(abs(v[i] - x1[i]) for v in verts) for i in 1:n]
    all(h .>= 0) || return nothing, "degenerate initial set"
    a = sqrt(n) .* max.(h, 1e-9)
    return LazySets.Ellipsoid(collect(float.(x1)), Matrix(LA.Diagonal(a .^ 2))), nothing
end

function forward_step!(
    provider,
    backend,
    opts::ForwardOptions,
    S,
    k::Int,
    E_k,
    xnext,
    uk,
    Qhat,
)
    @assert !isempty(opts.linearization_δu) "Set options.linearization_δu."

    xk = collect(LazySets.center(E_k))
    # The state box is known exactly from E_k — the whole point of forward mode.
    δx = opts.linearization_δx_margin .* _ellipsoid_axis_radii(E_k)
    δu = collect(Float64, opts.linearization_δu)

    approx = ST.build_affine_approximation(provider, xk, collect(uk); δx = δx, δu = δu)

    result = ST.solve_transition_forward(
        approx.system,
        E_k,
        collect(float.(xnext)),
        collect(uk),
        approx.Uformat,
        approx.Wformat,
        S,
        approx.lipschitz,
        backend;
        target_shape = opts.target_mode === :fixed ? Qhat : nothing,
        maxδu = opts.maxδu,
        λ = opts.λ,
        q_min = opts.q_min,
        q_max = opts.q_max,
        remainder_model = opts.remainder_model,
    )

    if !result.feasible
        return StepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            (; Xbar_radius = δx, Ubar_radius = δu, status = :forward_infeasible),
        )
    end

    # u-side consistency (the x-side holds by construction): the controller image
    # over E_k must stay inside the box the Hessian bound was taken on.
    uc, ru = _controller_image_axis_radii(result.controller, E_k)
    required_δu = abs.(uc .- collect(Float64, uk)) .+ ru
    if !all(required_δu .<= δu .+ 1e-8)
        return StepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            (;
                Xbar_radius = δx,
                Ubar_radius = δu,
                required_U_radius = required_δu,
                status = :forward_inconsistent_δu,
            ),
        )
    end

    E_next = result.target
    scale =
        Qhat === nothing ? NaN :
        LA.tr(Matrix{Float64}(LazySets.shape_matrix(E_next))) / LA.tr(Qhat)

    return StepRecord(
        k,
        :ok,
        result.cost,
        E_next,
        result.controller,
        (;
            Xbar_radius = δx,
            Ubar_radius = δu,
            required_U_radius = required_δu,
            contraction = scale,
            L = approx.lipschitz,
        ),
    )
end

"""
    ForwardCertifier(affine_provider, backend, options::ForwardOptions)

Ellipsoidal forward trajectory certifier: propagates a certified tube from the
entry ellipsoid along the nominal trajectory (see [`ForwardOptions`](@ref)).
Success requires every step certified, every enabled gate passed, and — when
`check_terminal` — the final tube inside the problem's target set. The result
reports the per-step contraction profile.
"""
mutable struct ForwardCertifier{AP, Backend, Opts} <: AbstractTrajectoryCertifier
    problem::Union{Nothing, PR.ProblemType}
    traj::Union{Nothing, ST.Trajectory}

    affine_provider::AP
    backend::Backend
    options::Opts

    result::Union{Nothing, CertificationResult}
    success::Bool
    solve_time_sec::Float64
end

function ForwardCertifier(affine_provider, backend, options::ForwardOptions)
    return ForwardCertifier(
        nothing,
        nothing,
        affine_provider,
        backend,
        options,
        nothing,
        false,
        0.0,
    )
end

function set_problem!(cert::ForwardCertifier, prob::PR.ProblemType)
    cert.problem = prob
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function set_trajectory!(cert::ForwardCertifier, traj::ST.Trajectory)
    cert.traj = traj
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

get_result(cert::ForwardCertifier) = cert.result
get_success(cert::ForwardCertifier) = cert.success
get_solve_time(cert::ForwardCertifier) = cert.solve_time_sec

function get_controller(cert::ForwardCertifier)
    res = cert.result
    (res === nothing || !res.success) && return nothing
    return ST.FunnelController(
        collect(res.lmi_data.kappas),
        collect(res.lmi_data.ellipsoids),
    )
end

function certify!(cert::ForwardCertifier)
    @assert cert.problem !== nothing "Call set_problem!(cert, problem) first."
    @assert cert.traj !== nothing "Call set_trajectory!(cert, trajectory) first."

    t0 = time()
    cert.result = _run_forward_chain(cert)
    cert.success = cert.result.success
    cert.solve_time_sec = time() - t0
    return cert
end

function _forward_gates(rec::StepRecord, problem, opts::ForwardOptions)
    rec.status == :ok || return rec
    reason = collapse_gate(rec, opts.r_min)
    reason === nothing || return _gate_failure(rec, :collapse, reason)
    if isfinite(opts.α_max) &&
       rec.summary.contraction isa Float64 &&
       !isnan(rec.summary.contraction) &&
       rec.summary.contraction > opts.α_max
        return _gate_failure(
            rec,
            :tube_inflation,
            "tube scale $(round(rec.summary.contraction; sigdigits = 3)) exceeds α_max",
        )
    end
    if opts.check_state_domain
        reason = state_domain_gate(rec, problem.system.X)
        reason === nothing || return _gate_failure(rec, :state_domain, reason)
    end
    return rec
end

function _run_forward_chain(cert::ForwardCertifier)
    t0 = time()
    problem = cert.problem
    opts = cert.options

    xs = collect(ST.states(cert.traj))
    us = collect(ST.inputs(cert.traj))
    K = length(us)
    @assert length(xs) == K + 1 "Expected length(xs) == length(us) + 1."

    nx = length(xs[1])
    nu = length(us[1])
    S = _transition_cost_matrix(opts.transition_cost, nx, nu)

    E_1, entry_reason = _entry_ellipsoid(problem, xs[1], opts)
    if E_1 === nothing
        return CertificationResult(
            false,
            0,
            Float64(time() - t0),
            StepRecord[],
            nothing,
            nothing,
            nothing,
            false,
            (; ellipsoids = [], kappas = [], entry_reason),
        )
    end
    Qhat =
        opts.target_mode === :fixed ? Matrix{Float64}(LazySets.shape_matrix(E_1)) : nothing

    domain_checked = opts.check_state_domain && _state_domain_supported(problem.system.X)

    steps = StepRecord[]
    ellipsoids = Any[E_1]
    E_k = E_1

    for k in 1:K
        rec = forward_step!(
            cert.affine_provider,
            cert.backend,
            opts,
            S,
            k,
            E_k,
            xs[k + 1],
            us[k],
            Qhat,
        )
        rec = _forward_gates(rec, problem, opts)
        push!(steps, rec)

        if rec.status != :ok
            kappas = _collect_kappas(steps)
            return CertificationResult(
                false,
                k,
                Float64(time() - t0),
                steps,
                nothing,
                nothing,
                1.0,   # entry coverage is exact by construction
                domain_checked,
                (; ellipsoids = copy(ellipsoids), kappas),
            )
        end

        E_k = rec.ellipsoid
        push!(ellipsoids, E_k)
    end

    terminal_contained =
        opts.check_terminal ? terminal_containment(E_k, problem.target_set) : nothing
    success = terminal_contained !== false
    kappas = _collect_kappas(steps)

    return CertificationResult(
        success,
        success ? nothing : K + 1,
        Float64(time() - t0),
        steps,
        success ? kappas : nothing,
        terminal_contained,
        1.0,
        domain_checked,
        (; ellipsoids = copy(ellipsoids), kappas),
    )
end
