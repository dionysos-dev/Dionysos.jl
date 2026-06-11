import LinearAlgebra as LA
using LaTeXStrings
import MathematicalSystems as MS
import Random
import StaticArrays: SVector

"""
    _eval_kappa_controller(kappa, x)

Evaluate a local controller `kappa` on the state `x`.
Supports affine maps, augmented matrices `[K ell]`, and the generic
`ST.get_control` fallback.
"""
function _eval_kappa_controller(kappa, x::AbstractVector)
    xvec = collect(Float64, x)

    if hasproperty(kappa, :A) && hasproperty(kappa, :b)
        A = Matrix{Float64}(getproperty(kappa, :A))
        b = vec(Float64.(getproperty(kappa, :b)))
        return collect(A * xvec + b)
    end
    if hasproperty(kappa, :A) && hasproperty(kappa, :c)
        A = Matrix{Float64}(getproperty(kappa, :A))
        b = vec(Float64.(getproperty(kappa, :c)))
        return collect(A * xvec + b)
    end
    if kappa isa AbstractMatrix
        K = Matrix{Float64}(kappa[:, 1:(end - 1)])
        ell = vec(Float64.(kappa[:, end]))
        return collect(K * xvec + ell)
    end

    return collect(Float64.(ST.get_control(kappa, xvec)))
end

function _sample_points_uniform_in_hyperrectangle(
    rect::UT.HyperRectangle,
    n::Int;
    rng = Random.default_rng(),
)
    n >= 1 || error("n must be >= 1.")

    lb = collect(Float64, rect.lb)
    ub = collect(Float64, rect.ub)
    nx = length(lb)
    pts = Vector{Vector{Float64}}(undef, n)

    for i in 1:n
        pts[i] = [lb[j] + Random.rand(rng) * (ub[j] - lb[j]) for j in 1:nx]
    end

    return pts
end

function _sample_points_uniform_in_ellipsoid(
    E::UT.Ellipsoid,
    n::Int;
    rng = Random.default_rng(),
    tol::Real = 1.0e-8,
)
    n >= 1 || error("n must be >= 1.")

    nx = length(E.c)
    U = LA.cholesky(LA.Symmetric(Matrix{Float64}(E.P))).U
    c = vec(Float64.(E.c))
    pts = Vector{Vector{Float64}}(undef, n)

    for i in 1:n
        v = Random.randn(rng, nx)
        nv = LA.norm(v)
        while nv <= 1.0e-14
            v = Random.randn(rng, nx)
            nv = LA.norm(v)
        end
        rho = Random.rand(rng)^(1 / nx)
        z = (rho / nv) .* v
        pts[i] = collect(c .+ U \ z)
    end

    validation = _validate_ellipsoid_samples(E, pts; tol = tol)
    @assert validation.max_rho <= 1.0 + tol
    return pts
end

function _ellipsoid_rho(E::UT.Ellipsoid, x::AbstractVector)
    dx = collect(Float64, x) .- vec(Float64.(E.c))
    return Float64(dx' * Matrix{Float64}(E.P) * dx)
end

"""
    inflated_ellipsoid(E, alpha)

Return the ellipsoid obtained by multiplying the radii of `E` by `alpha`.
Dionysos stores ellipsoids as `(x-c)'P(x-c) <= 1`, so radius inflation by
`alpha` is represented by the shape matrix `P / alpha^2`.
"""
function inflated_ellipsoid(E::UT.Ellipsoid, alpha::Real)
    isfinite(alpha) && alpha > 0.0 || error("alpha must be positive and finite.")
    return UT.Ellipsoid(Matrix{Float64}(E.P) ./ (Float64(alpha)^2), vec(Float64.(E.c)))
end

"""
    sample_points_uniform_in_inflated_ellipsoid(E, n; alpha, rng)

Sample uniformly in the ellipsoid whose radii are `alpha` times the radii of
`E`.
"""
function sample_points_uniform_in_inflated_ellipsoid(
    E::UT.Ellipsoid,
    n::Int;
    alpha::Real = 1.0,
    rng = Random.default_rng(),
)
    return _sample_points_uniform_in_ellipsoid(inflated_ellipsoid(E, alpha), n; rng = rng)
end

"""
    sample_points_uniform_in_ellipsoid_shell(E, n; alpha_inner, alpha_outer, rng)

Sample uniformly in the shell between radius inflations `alpha_inner` and
`alpha_outer` of `E`.  The shell is useful for empirical stress tests because
it avoids statistics being dominated by points close to the ellipsoid center.
"""
function sample_points_uniform_in_ellipsoid_shell(
    E::UT.Ellipsoid,
    n::Int;
    alpha_inner::Real,
    alpha_outer::Real,
    rng = Random.default_rng(),
)
    n >= 1 || error("n must be >= 1.")
    isfinite(alpha_inner) && alpha_inner >= 0.0 ||
        error("alpha_inner must be nonnegative and finite.")
    isfinite(alpha_outer) && alpha_outer > alpha_inner ||
        error("alpha_outer must be finite and larger than alpha_inner.")

    nx = length(E.c)
    U = LA.cholesky(LA.Symmetric(Matrix{Float64}(E.P))).U
    c = vec(Float64.(E.c))
    pts = Vector{Vector{Float64}}(undef, n)
    inner_pow = Float64(alpha_inner)^nx
    outer_pow = Float64(alpha_outer)^nx

    for i in 1:n
        v = Random.randn(rng, nx)
        nv = LA.norm(v)
        while nv <= 1.0e-14
            v = Random.randn(rng, nx)
            nv = LA.norm(v)
        end
        rho = (inner_pow + Random.rand(rng) * (outer_pow - inner_pow))^(1 / nx)
        z = (rho / nv) .* v
        pts[i] = collect(c .+ U \ z)
    end

    outer_E = inflated_ellipsoid(E, alpha_outer)
    validation = _validate_ellipsoid_samples(outer_E, pts)
    @assert validation.n_outside_ellipsoid == 0
    return pts
end

function _validate_ellipsoid_samples(
    E::UT.Ellipsoid,
    pts;
    tol::Real = 1.0e-8,
    initial_set = nothing,
)
    rhos = [_ellipsoid_rho(E, x) for x in pts]
    max_rho = isempty(rhos) ? 0.0 : maximum(rhos)
    max_violation = max(0.0, max_rho - 1.0)
    n_outside_ellipsoid = count(rho -> rho > 1.0 + tol, rhos)
    n_outside_initial_set =
        initial_set === nothing ? nothing : count(x -> !(x ∈ initial_set), pts)
    return (;
        n_samples = length(pts),
        max_rho = max_rho,
        max_ellipsoid_violation = max_violation,
        n_outside_ellipsoid = n_outside_ellipsoid,
        n_outside_initial_set = n_outside_initial_set,
    )
end

function _validate_initial_samples(
    sample_set,
    pts;
    initial_set = nothing,
    tol::Real = 1.0e-8,
)
    if sample_set isa UT.Ellipsoid
        return _validate_ellipsoid_samples(
            sample_set,
            pts;
            tol = tol,
            initial_set = initial_set,
        )
    end

    n_outside_initial_set =
        initial_set === nothing ? nothing : count(x -> !(x ∈ initial_set), pts)
    return (;
        n_samples = length(pts),
        max_rho = nothing,
        max_ellipsoid_violation = nothing,
        n_outside_ellipsoid = nothing,
        n_outside_initial_set = n_outside_initial_set,
    )
end

"""
    sample_points_uniform_in_set(set, n; rng)

Sample `n` points uniformly in a supported set.
Currently supports `UT.HyperRectangle` and `UT.Ellipsoid`.
"""
function sample_points_uniform_in_set(set, n::Int; rng = Random.default_rng())
    set isa UT.HyperRectangle &&
        return _sample_points_uniform_in_hyperrectangle(set, n; rng = rng)
    set isa UT.Ellipsoid && return _sample_points_uniform_in_ellipsoid(set, n; rng = rng)

    return error("Unsupported sampling set type: $(typeof(set)).")
end

function _project_input_to_domain(u::AbstractVector, u_domain)
    u_domain isa UT.HyperRectangle || return collect(Float64, u)
    return collect(clamp.(u, u_domain.lb, u_domain.ub))
end

function _input_in_domain(u::AbstractVector, u_domain)
    try
        return u ∈ u_domain
    catch
        return true
    end
end

function _build_certified_kappa_chain(cert_result)
    cert_result === nothing && error("No certification result available.")
    hasproperty(cert_result, :steps) || error("Certification result has no `steps` field.")

    valid_steps = [
        step for step in cert_result.steps if
        step.status == :ok && step.kappa !== nothing && step.ellipsoid !== nothing
    ]
    isempty(valid_steps) &&
        error("No certified kappa controller found in the certification result.")

    sort!(valid_steps; by = step -> step.k)
    k_sequence = [step.k for step in valid_steps]
    expected_k = collect(first(k_sequence):last(k_sequence))
    k_sequence == expected_k || error(
        "The kappa chain is not contiguous: got $(k_sequence), expected $(expected_k).",
    )
    first(k_sequence) == 1 || error(
        "The kappa chain does not start at k = 1, so it cannot be replayed from the initial set.",
    )

    hasproperty(cert_result, :lmi_data) &&
    cert_result.lmi_data !== nothing &&
    hasproperty(cert_result.lmi_data, :ellipsoids) &&
    !isempty(cert_result.lmi_data.ellipsoids) || error(
        "The certification result does not expose its ellipsoid chain in `lmi_data.ellipsoids`.",
    )

    kappa_by_k = Dict{Int, Any}()
    ellipsoid_by_k = Dict{Int, UT.Ellipsoid}()
    for step in valid_steps
        kappa_by_k[step.k] = step.kappa
        ellipsoid_by_k[step.k] = step.ellipsoid
    end

    # `run_backward_chain!` stores ellipsoids as [E_terminal, ..., E_initial].
    # The step records already contain E_k.  The only missing element for
    # replay is the terminal ellipsoid E_{K+1}.
    k_terminal = last(k_sequence) + 1
    ellipsoid_by_k[k_terminal] = first(cert_result.lmi_data.ellipsoids)

    return (;
        k_sequence = k_sequence,
        kappa_by_k = kappa_by_k,
        ellipsoid_by_k = ellipsoid_by_k,
        initial_ellipsoid = ellipsoid_by_k[1],
        target_ellipsoid = ellipsoid_by_k[k_terminal],
    )
end

function _domain_set(system)
    !hasproperty(system, :X) && return nothing
    X = system.X
    return X isa UT.LazySetMinus ? X.A : X
end

function _obstacle_set(system)
    if hasproperty(system, :obstacles)
        obs = system.obstacles
        isempty(obs) && return nothing
        return obs
    end

    if hasproperty(system, :X) && system.X isa UT.LazySetMinus
        return system.X.B
    end

    return nothing
end

function _point_in_single_set(x, set)
    try
        return x ∈ set
    catch
        if hasproperty(set, :sets)
            return any(s -> _point_in_single_set(x, s), set.sets)
        end
        if hasproperty(set, :lb)
            d = length(set.lb)
            d <= length(x) && return x[1:d] ∈ set
        end
        rethrow()
    end
end

function _point_in_any_set(x, set)
    set === nothing && return false
    if set isa AbstractVector
        return any(s -> _point_in_single_set(x, s), set)
    end
    return _point_in_single_set(x, set)
end

function _build_periodic_maps(;
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
)
    if periodic_dims === nothing || periodic_periods === nothing
        wrap_state = x -> collect(Float64, x)
        lift_state_near_reference = (x, _) -> collect(Float64, x)
        return (; has_periodic = false, wrap_state, lift_state_near_reference)
    end

    start =
        periodic_start === nothing ? zeros(length(periodic_dims)) :
        collect(Float64, periodic_start)

    wrap_state = function (x)
        xv = collect(Float64, x)
        wrap_vector_periodic!(xv, periodic_dims, periodic_periods, start)
        return xv
    end

    lift_state_near_reference = function (x, xref)
        xv = collect(Float64, x)
        y = copy(xv)

        for i in eachindex(periodic_dims)
            d = Int(periodic_dims[i])
            p = Float64(periodic_periods[i])
            y[d] += round((Float64(xref[d]) - y[d]) / p) * p
        end

        return y
    end

    return (; has_periodic = true, wrap_state, lift_state_near_reference)
end

function _build_discrete_system_map(system, Ts; num_substeps::Int = 1)
    num_substeps >= 1 || error("num_substeps must be >= 1.")

    try
        disc = ST.discretize_continuous_system(
            system,
            Float64(Ts);
            num_substeps = num_substeps,
        )
        return MS.mapping(disc)
    catch
        return MS.mapping(system)
    end
end

function _rollout_label(stats)
    stats.certified_success && return "certified_success"
    stats.closed_loop_success && return "target_reached"
    stats.hit_obstacle && return "obstacle"
    stats.left_domain && return "domain_exit"
    stats.input_violation && return "input_violation"
    stats.first_exit_k !== nothing && return "ellipsoid_exit"
    stats.started_outside_initial_ellipsoid && return "starts_outside"
    stats.final_in_target || return "miss_target"
    return "other"
end

function _rollout_color(stats)
    return stats.closed_loop_success ? PLOT_COLORS[:rollout] : PLOT_COLORS[:failure]
end

_rollout_is_success(stats) = stats.closed_loop_success

_rollout_line_alpha(stats, default_alpha::Real) =
    _rollout_is_success(stats) ? min(Float64(default_alpha), 0.20) : 0.88

_rollout_linewidth(stats, default_lw::Real) =
    _rollout_is_success(stats) ? Float64(default_lw) : max(Float64(default_lw), 1.15)

function _rollout_pretty_label(lbl::String)
    lbl == "certified_success" && return "Certified success"
    lbl == "target_reached" && return "Empirical success"
    lbl == "obstacle" && return "Obstacle hit"
    lbl == "domain_exit" && return "Domain exit"
    lbl == "input_violation" && return "Input violation"
    lbl == "ellipsoid_exit" && return "Certificate exit"
    lbl == "starts_outside" && return "Starts outside initial ellipsoid"
    lbl == "miss_target" && return "Target missed"
    return "Other"
end

_show_rollout_label_in_legend(lbl::String) = !(lbl in ("certified_success", "domain_exit"))

function _resolve_periodic_plot_config(
    stat_result;
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
)
    cfg =
        hasproperty(stat_result, :config) ? getproperty(stat_result, :config) :
        (
            hasproperty(stat_result, :run_config) ? getproperty(stat_result, :run_config) :
            nothing
        )

    if periodic_dims === nothing && cfg !== nothing && hasproperty(cfg, :periodic_dims)
        periodic_dims = cfg.periodic_dims
    end
    if periodic_periods === nothing &&
       cfg !== nothing &&
       hasproperty(cfg, :periodic_periods)
        periodic_periods = cfg.periodic_periods
    end
    if periodic_start === nothing && cfg !== nothing && hasproperty(cfg, :periodic_start)
        periodic_start = cfg.periodic_start
    end

    return (; periodic_dims, periodic_periods, periodic_start)
end

function _transform_state_list_for_plot(
    state_list;
    unwrap_angles::Bool,
    wrap_angles::Bool,
    periodic_dims,
    periodic_periods,
    periodic_start,
)
    xs = [SVector{length(x), Float64}(x) for x in state_list]

    if unwrap_angles && periodic_dims !== nothing && periodic_periods !== nothing
        xs = unwrap_periodic_state_list(xs, periodic_dims, periodic_periods)
    end
    if wrap_angles &&
       periodic_dims !== nothing &&
       periodic_periods !== nothing &&
       periodic_start !== nothing
        xs = wrap_periodic_state_list(xs, periodic_dims, periodic_periods, periodic_start)
    end

    return xs
end

function _transform_single_state_for_plot(
    x;
    unwrap_angles::Bool,
    wrap_angles::Bool,
    periodic_dims,
    periodic_periods,
    periodic_start,
)
    xs = _transform_state_list_for_plot(
        [x];
        unwrap_angles = unwrap_angles,
        wrap_angles = wrap_angles,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    return xs[1]
end

function _ellipsoids_for_plot(stat_result; max_keep::Int = 60)
    ks = sort(collect(keys(stat_result.ellipsoid_by_k)))
    ells = [stat_result.ellipsoid_by_k[k] for k in ks]

    length(ells) <= max_keep && return ells

    stride = max(1, cld(length(ells), max_keep))
    subset = ells[1:stride:end]
    last(subset) === last(ells) && return subset
    push!(subset, last(ells))
    return subset
end

function _ellipsoid_records_for_plot(stat_result; max_keep::Int = 60)
    ks = sort(collect(keys(stat_result.ellipsoid_by_k)))
    records = [(; state_index = k, ellipsoid = stat_result.ellipsoid_by_k[k]) for k in ks]

    length(records) <= max_keep && return records

    stride = max(1, cld(length(records), max_keep))
    subset = records[1:stride:end]
    last(subset).state_index == last(records).state_index && return subset
    push!(subset, last(records))
    return subset
end

"""
    run_kappa_statistical_check(run_result; kwargs...)

Sample points in the initial set, replay the certified kappa controllers on
the concrete system, track obstacle/domain violations, target reachability,
and the first ellipsoid exit for each rollout.

Important:
- this function expects a successful ellipsoidal certification,
- points are sampled from the first certified ellipsoid by default,
- for periodic coordinates, pass `periodic_dims`, `periodic_periods`,
  and `periodic_start` so that the state is wrapped for the concrete
  simulation and lifted near the certification trajectory for ellipsoid
  and kappa evaluation.
"""
function run_kappa_statistical_check(
    run_result;
    n_samples::Int = 100,
    sample_set = nothing,
    sample_from::Symbol = :initial_ellipsoid,
    x0_samples = nothing,
    num_substeps::Int = 1,
    project_inputs::Bool = false,
    check_domain::Bool = true,
    domain_set = nothing,
    obstacle_set = nothing,
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
    stop_on_target::Bool = true,
    verbose::Bool = true,
    rng = Random.default_rng(),
)
    if x0_samples !== nothing
        n_samples = length(x0_samples)
    end
    n_samples >= 1 || error("n_samples must be >= 1.")
    num_substeps >= 1 || error("num_substeps must be >= 1.")

    hasproperty(run_result, :problem) ||
        error("`run_result` must expose a `problem` field.")
    hasproperty(run_result, :result) || error("`run_result` must expose a `result` field.")
    hasproperty(run_result, :certification_candidate) ||
        error("`run_result` must expose a `certification_candidate` field.")

    problem = run_result.problem
    cert_result = run_result.result.certification
    chain = _build_certified_kappa_chain(cert_result)

    cert_candidate = run_result.certification_candidate
    cert_xs = collect(ST.enum_elems(cert_candidate.x_traj))
    length(cert_xs) >= last(chain.k_sequence) + 1 ||
        error("The certification trajectory is shorter than the certified kappa chain.")

    if hasproperty(run_result, :config)
        cfg = run_result.config
        if periodic_dims === nothing && hasproperty(cfg, :periodic_dims)
            periodic_dims = cfg.periodic_dims
        end
        if periodic_periods === nothing && hasproperty(cfg, :periodic_periods)
            periodic_periods = cfg.periodic_periods
        end
        if periodic_start === nothing && hasproperty(cfg, :periodic_start)
            periodic_start = cfg.periodic_start
        end
    end

    if sample_set === nothing
        if sample_from == :initial_ellipsoid
            sample_set = chain.initial_ellipsoid
        elseif sample_from == :initial_set
            sample_set = problem.initial_set
        else
            error("Unsupported `sample_from` mode: $(sample_from).")
        end
    end

    if domain_set === nothing
        domain_set = _domain_set(problem.system)
    end
    if obstacle_set === nothing
        obstacle_set = _obstacle_set(problem.system)
    end

    maps = _build_periodic_maps(;
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    target_reference =
        hasproperty(problem, :target_set) ?
        collect(Float64, UT.get_center(problem.target_set)) : nothing
    target_contains = function (x)
        target_reference === nothing && return false
        x_phys = maps.wrap_state(x)
        x_target = maps.lift_state_near_reference(x_phys, target_reference)
        return x_target ∈ problem.target_set
    end
    f_disc = _build_discrete_system_map(
        problem.system,
        cert_candidate.Ts;
        num_substeps = num_substeps,
    )
    x0_samples =
        x0_samples === nothing ?
        sample_points_uniform_in_set(sample_set, n_samples; rng = rng) :
        [collect(Float64, x0) for x0 in x0_samples]
    sampling_validation =
        _validate_initial_samples(sample_set, x0_samples; initial_set = problem.initial_set)
    if sampling_validation.n_outside_ellipsoid !== nothing
        @assert sampling_validation.n_outside_ellipsoid == 0
    end

    x_rollouts = Vector{Vector{Vector{Float64}}}(undef, n_samples)
    u_rollouts = Vector{Vector{Vector{Float64}}}(undef, n_samples)
    rollout_stats = Vector{NamedTuple}(undef, n_samples)
    ellipsoid_exit_points = NamedTuple[]

    for s in 1:n_samples
        x = maps.wrap_state(x0_samples[s])
        x_hist = Vector{Vector{Float64}}([copy(x)])
        u_hist = Vector{Vector{Float64}}()

        hit_obstacle = false
        left_domain = false
        input_violation = false
        first_input_violation_k = nothing
        target_reached = false
        target_reach_k = nothing
        final_in_target = false
        first_exit_k = nothing
        first_exit_phase = nothing
        first_exit_state = nothing
        x_init_chain = maps.lift_state_near_reference(x, cert_xs[first(chain.k_sequence)])
        inside_initial_ellipsoid = x_init_chain ∈ chain.initial_ellipsoid
        started_outside_initial_ellipsoid = !inside_initial_ellipsoid
        ever_inside_ellipsoid_chain = inside_initial_ellipsoid

        for k in chain.k_sequence
            x_phys = maps.wrap_state(x)
            x_chain = maps.lift_state_near_reference(x_phys, cert_xs[k])
            inside_current = x_chain ∈ chain.ellipsoid_by_k[k]

            if !ever_inside_ellipsoid_chain && inside_current
                ever_inside_ellipsoid_chain = true
            end

            if ever_inside_ellipsoid_chain && !inside_current && first_exit_k === nothing
                first_exit_k = k
                first_exit_phase = :before_step
                first_exit_state = copy(x_phys)
            end

            u = _eval_kappa_controller(chain.kappa_by_k[k], x_chain)
            if hasproperty(problem.system, :U) && !_input_in_domain(u, problem.system.U)
                input_violation = true
                first_input_violation_k === nothing && (first_input_violation_k = k)
                if !project_inputs
                    push!(u_hist, copy(u))
                    break
                end
            end
            if project_inputs && hasproperty(problem.system, :U)
                u = _project_input_to_domain(u, problem.system.U)
            end
            push!(u_hist, copy(u))

            x_next = maps.wrap_state(f_disc(x_phys, u))
            push!(x_hist, copy(x_next))

            x_next_chain = maps.lift_state_near_reference(x_next, cert_xs[k + 1])
            inside_next = x_next_chain ∈ chain.ellipsoid_by_k[k + 1]

            if !ever_inside_ellipsoid_chain && inside_next
                ever_inside_ellipsoid_chain = true
            end

            if ever_inside_ellipsoid_chain && !inside_next && first_exit_k === nothing
                first_exit_k = k + 1
                first_exit_phase = :after_step
                first_exit_state = copy(x_next)
            end

            if _point_in_any_set(x_next, obstacle_set)
                hit_obstacle = true
                x = x_next
                break
            end

            if check_domain && domain_set !== nothing && !(x_next ∈ domain_set)
                left_domain = true
                x = x_next
                break
            end

            if target_contains(x_next)
                target_reached = true
                target_reach_k = k + 1
                x = x_next
                stop_on_target && break
            end

            x = x_next
        end

        final_state = maps.wrap_state(x)
        final_in_target = target_contains(final_state)
        target_reached |= final_in_target
        closed_loop_success =
            !hit_obstacle && !left_domain && !input_violation && target_reached
        certified_success =
            !hit_obstacle &&
            !left_domain &&
            !input_violation &&
            inside_initial_ellipsoid &&
            first_exit_k === nothing &&
            target_reached

        x_rollouts[s] = x_hist
        u_rollouts[s] = u_hist
        rollout_stats[s] = (;
            success = certified_success,
            certified_success = certified_success,
            closed_loop_success = closed_loop_success,
            hit_obstacle = hit_obstacle,
            left_domain = left_domain,
            input_violation = input_violation,
            first_input_violation_k = first_input_violation_k,
            target_reached = target_reached,
            target_reach_k = target_reach_k,
            final_in_target = final_in_target,
            inside_initial_ellipsoid = inside_initial_ellipsoid,
            started_outside_initial_ellipsoid = started_outside_initial_ellipsoid,
            first_exit_k = first_exit_k,
            first_exit_phase = first_exit_phase,
            first_exit_state = first_exit_state,
        )

        if first_exit_k !== nothing
            push!(
                ellipsoid_exit_points,
                (;
                    sample_idx = s,
                    k = first_exit_k,
                    phase = first_exit_phase,
                    state = first_exit_state,
                ),
            )
        end
    end

    n_certified_success = count(r -> r.certified_success, rollout_stats)
    n_closed_loop_success = count(r -> r.closed_loop_success, rollout_stats)
    n_obstacle = count(r -> r.hit_obstacle, rollout_stats)
    n_left_domain = count(r -> r.left_domain, rollout_stats)
    n_input_violation = count(r -> r.input_violation, rollout_stats)
    n_target = count(r -> r.final_in_target, rollout_stats)
    n_ellipsoid_exit = count(r -> r.first_exit_k !== nothing, rollout_stats)
    n_exit_before = count(r -> r.first_exit_phase == :before_step, rollout_stats)
    n_exit_after = count(r -> r.first_exit_phase == :after_step, rollout_stats)
    n_started_outside_initial_ellipsoid =
        count(r -> r.started_outside_initial_ellipsoid, rollout_stats)

    exit_histogram = Dict{Int, Int}()
    for stats in rollout_stats
        stats.first_exit_k === nothing && continue
        exit_histogram[stats.first_exit_k] = get(exit_histogram, stats.first_exit_k, 0) + 1
    end

    summary = (;
        n_samples = n_samples,
        n_success = n_certified_success,
        n_certified_success = n_certified_success,
        n_closed_loop_success = n_closed_loop_success,
        n_hit_obstacle = n_obstacle,
        n_left_domain = n_left_domain,
        n_input_violation = n_input_violation,
        n_final_in_target = n_target,
        n_ellipsoid_exit = n_ellipsoid_exit,
        n_exit_before_step = n_exit_before,
        n_exit_after_step = n_exit_after,
        n_started_outside_initial_ellipsoid = n_started_outside_initial_ellipsoid,
        n_outside_initial_ellipsoid = n_started_outside_initial_ellipsoid,
        sampling_max_initial_ellipsoid_violation = sampling_validation.max_ellipsoid_violation,
        sampling_n_outside_initial_ellipsoid = sampling_validation.n_outside_ellipsoid,
        sampling_n_outside_initial_set = sampling_validation.n_outside_initial_set,
        success_rate = n_certified_success / n_samples,
        certified_success_rate = n_certified_success / n_samples,
        closed_loop_success_rate = n_closed_loop_success / n_samples,
        obstacle_rate = n_obstacle / n_samples,
        domain_exit_rate = n_left_domain / n_samples,
        input_violation_rate = n_input_violation / n_samples,
        final_target_rate = n_target / n_samples,
        ellipsoid_exit_rate = n_ellipsoid_exit / n_samples,
        initial_ellipsoid_coverage_rate = 1.0 -
                                          n_started_outside_initial_ellipsoid / n_samples,
    )

    stat_result = (;
        problem = problem,
        run_config = hasproperty(run_result, :config) ? run_result.config : nothing,
        config = hasproperty(run_result, :config) ? run_result.config : nothing,
        outputs = hasproperty(run_result, :outputs) ? run_result.outputs : nothing,
        sample_set = sample_set,
        sample_from = sample_from,
        domain_set = domain_set,
        obstacle_set = obstacle_set,
        certification_result = cert_result,
        certification_candidate = cert_candidate,
        k_sequence = chain.k_sequence,
        kappa_by_k = chain.kappa_by_k,
        ellipsoid_by_k = chain.ellipsoid_by_k,
        initial_ellipsoid = chain.initial_ellipsoid,
        target_ellipsoid = chain.target_ellipsoid,
        reference_states = cert_xs,
        x0_samples = x0_samples,
        sampling_validation = sampling_validation,
        x_rollouts = x_rollouts,
        u_rollouts = u_rollouts,
        rollout_stats = rollout_stats,
        ellipsoid_exit_points = ellipsoid_exit_points,
        exit_histogram = exit_histogram,
        summary = summary,
    )

    verbose && print_kappa_statistical_summary(stat_result)
    return stat_result
end

"""
    run_closed_loop_rollouts_from_initial_states(run_result, x0_samples; kwargs...)

Replay the certified feedback chain from an explicit batch of initial states.
This is a thin wrapper around `run_kappa_statistical_check`; it is intended for
empirical diagnostics where the sampling distribution is chosen outside the
standard certified initial ellipsoid.
"""
function run_closed_loop_rollouts_from_initial_states(
    run_result,
    x0_samples;
    sample_set = nothing,
    kwargs...,
)
    isempty(x0_samples) && error("x0_samples must not be empty.")
    return run_kappa_statistical_check(
        run_result;
        n_samples = length(x0_samples),
        sample_set = sample_set,
        x0_samples = x0_samples,
        kwargs...,
    )
end

function _empirical_inflation_summary_row(stat_result, alpha_inner::Real, alpha_outer::Real)
    stats = stat_result.rollout_stats
    n = length(stats)
    n >= 1 || error("A stress-test row needs at least one rollout.")

    success_count = count(r -> r.closed_loop_success, stats)
    target_reached_count = count(r -> r.target_reached, stats)
    obstacle_violation_count = count(r -> r.hit_obstacle, stats)
    state_violation_count = count(r -> r.left_domain, stats)
    input_violation_count = count(r -> r.input_violation, stats)
    chain_violation_count =
        count(r -> r.started_outside_initial_ellipsoid || r.first_exit_k !== nothing, stats)
    started_outside_initial_ellipsoid_count =
        count(r -> r.started_outside_initial_ellipsoid, stats)

    return (;
        alpha_inner = Float64(alpha_inner),
        alpha_outer = Float64(alpha_outer),
        n_samples = n,
        success_count = success_count,
        success_rate = success_count / n,
        target_reached_count = target_reached_count,
        target_reached_rate = target_reached_count / n,
        obstacle_violation_count = obstacle_violation_count,
        obstacle_violation_rate = obstacle_violation_count / n,
        state_violation_count = state_violation_count,
        state_violation_rate = state_violation_count / n,
        input_violation_count = input_violation_count,
        input_violation_rate = input_violation_count / n,
        chain_violation_count = chain_violation_count,
        chain_violation_rate = chain_violation_count / n,
        started_outside_initial_ellipsoid_count = started_outside_initial_ellipsoid_count,
        started_outside_initial_ellipsoid_rate = started_outside_initial_ellipsoid_count /
                                                 n,
    )
end

function print_empirical_inflation_stress_summary(summary_rows)
    println("\n=== Empirical inflation stress test ===")
    println("These statistics are empirical only; outside E0 they are not certified.")
    for row in summary_rows
        label =
            row.alpha_inner == 0.0 ? "E0($(row.alpha_outer))" :
            "shell $(row.alpha_inner) -> $(row.alpha_outer)"
        println(
            label,
            ": success = ",
            row.success_count,
            "/",
            row.n_samples,
            " (",
            round(100 * row.success_rate; digits = 1),
            "%), obstacle = ",
            round(100 * row.obstacle_violation_rate; digits = 1),
            "%, state = ",
            round(100 * row.state_violation_rate; digits = 1),
            "%, input = ",
            round(100 * row.input_violation_rate; digits = 1),
            "%, chain exit = ",
            round(100 * row.chain_violation_rate; digits = 1),
            "%",
        )
    end
    return nothing
end

"""
    run_empirical_inflation_stress_test(run_result; kwargs...)

Empirically tests the certified feedback sequence on inflated versions of the
certified initial ellipsoid.  At `alpha = 1`, the samples are in `E0`.  For
`alpha > 1`, and especially in shell mode, the samples are outside the formal
guarantee; the result is a basin-of-success diagnostic, not a certificate.
"""
function run_empirical_inflation_stress_test(
    run_result;
    inflation_factors = [1.0, 1.1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0],
    samples_per_alpha::Int = 2000,
    shell_sampling::Bool = true,
    selected_rollout_alpha::Union{Nothing, Real} = 2.0,
    num_substeps::Int = 1,
    project_inputs::Bool = false,
    check_domain::Bool = true,
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
    stop_on_target::Bool = true,
    verbose::Bool = true,
    rng = Random.default_rng(),
)
    samples_per_alpha >= 1 || error("samples_per_alpha must be >= 1.")
    alphas = sort(unique(Float64.(inflation_factors)))
    isempty(alphas) && error("inflation_factors must not be empty.")
    all(alpha -> isfinite(alpha) && alpha >= 1.0, alphas) ||
        error("inflation_factors must be finite and >= 1.")
    first(alphas) == 1.0 || pushfirst!(alphas, 1.0)

    cert_result = run_result.result.certification
    chain = _build_certified_kappa_chain(cert_result)
    E0 = chain.initial_ellipsoid

    summary_rows = NamedTuple[]
    rollout_results = NamedTuple[]
    selected_rollout_result = nothing

    for (idx, alpha_outer) in enumerate(alphas)
        alpha_inner = shell_sampling && idx > 1 ? alphas[idx - 1] : 0.0
        x0_samples =
            shell_sampling && idx > 1 ?
            sample_points_uniform_in_ellipsoid_shell(
                E0,
                samples_per_alpha;
                alpha_inner = alpha_inner,
                alpha_outer = alpha_outer,
                rng = rng,
            ) :
            sample_points_uniform_in_inflated_ellipsoid(
                E0,
                samples_per_alpha;
                alpha = alpha_outer,
                rng = rng,
            )
        outer_set = inflated_ellipsoid(E0, alpha_outer)

        stat_result = run_closed_loop_rollouts_from_initial_states(
            run_result,
            x0_samples;
            sample_set = outer_set,
            sample_from = :empirical_inflation_stress,
            num_substeps = num_substeps,
            project_inputs = project_inputs,
            check_domain = check_domain,
            periodic_dims = periodic_dims,
            periodic_periods = periodic_periods,
            periodic_start = periodic_start,
            stop_on_target = stop_on_target,
            verbose = false,
        )
        row = _empirical_inflation_summary_row(stat_result, alpha_inner, alpha_outer)
        push!(summary_rows, row)
        push!(rollout_results, (; alpha_inner, alpha_outer, stat_result))

        if selected_rollout_alpha !== nothing && isapprox(
            alpha_outer,
            Float64(selected_rollout_alpha);
            atol = 1.0e-10,
            rtol = 0.0,
        )
            selected_rollout_result = stat_result
        end
    end

    verbose && print_empirical_inflation_stress_summary(summary_rows)
    return (;
        summary_rows = summary_rows,
        rollout_results = rollout_results,
        selected_rollout_result = selected_rollout_result,
        initial_ellipsoid = E0,
        inflation_factors = alphas,
        shell_sampling = shell_sampling,
        samples_per_alpha = samples_per_alpha,
        note = "Empirical stress test only; alpha > 1 is outside the certified guarantee.",
    )
end

function _csv_value(x)
    x === nothing && return ""
    x isa AbstractFloat && return string(x)
    return string(x)
end

function save_empirical_inflation_stress_csv(path::AbstractString, summary_rows)
    mkpath(dirname(path))
    columns = [
        :alpha_inner,
        :alpha_outer,
        :n_samples,
        :success_count,
        :success_rate,
        :target_reached_count,
        :target_reached_rate,
        :obstacle_violation_count,
        :obstacle_violation_rate,
        :state_violation_count,
        :state_violation_rate,
        :input_violation_count,
        :input_violation_rate,
        :chain_violation_count,
        :chain_violation_rate,
        :started_outside_initial_ellipsoid_count,
        :started_outside_initial_ellipsoid_rate,
    ]

    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for row in summary_rows
            println(io, join((_csv_value(getproperty(row, col)) for col in columns), ","))
        end
    end
    return path
end

function plot_empirical_inflation_stress_rates(
    summary_rows;
    output_dir::AbstractString = pwd(),
    filename::AbstractString = "empirical_success_vs_inflation.pdf",
)
    isempty(summary_rows) && error("No empirical stress-test rows to plot.")
    mkpath(output_dir)

    alpha = [row.alpha_outer for row in summary_rows]
    success = [row.success_rate for row in summary_rows]
    obstacle = [row.obstacle_violation_rate for row in summary_rows]
    state = [row.state_violation_rate for row in summary_rows]
    input = [row.input_violation_rate for row in summary_rows]

    fig = plot(
        alpha,
        success;
        thesis_plot_kwargs(; legend = :right, size = (760, 500))...,
        xlabel = L"\alpha",
        ylabel = "empirical rate",
        ylims = (0.0, 1.0),
        color = PLOT_COLORS[:nominal],
        lw = 2.6,
        marker = :circle,
        ms = 4.0,
        label = "success",
    )
    plot!(
        fig,
        alpha,
        obstacle;
        color = PLOT_COLORS[:failure],
        lw = 2.0,
        marker = :diamond,
        ms = 3.8,
        label = "obstacle violation",
    )
    plot!(
        fig,
        alpha,
        state;
        color = PLOT_COLORS[:constraint],
        lw = 2.0,
        ls = :dash,
        marker = :utriangle,
        ms = 3.8,
        label = "state-domain violation",
    )
    plot!(
        fig,
        alpha,
        input;
        color = PLOT_COLORS[:state2],
        lw = 2.0,
        ls = :dot,
        marker = :rect,
        ms = 3.5,
        label = "input-domain violation",
    )

    out_path = joinpath(output_dir, filename)
    savefig(fig, out_path)
    display(fig)
    return out_path
end

"""
    print_kappa_statistical_summary(stat_result)

Print a compact summary of the empirical kappa replay.
"""
function print_kappa_statistical_summary(stat_result)
    summary = stat_result.summary

    println("\n=== Empirical closed-loop validation ===")
    println("samples                  = ", summary.n_samples)
    println(
        "certified success        = ",
        summary.n_certified_success,
        " (",
        round(100 * summary.certified_success_rate; digits = 1),
        "%)",
    )
    println(
        "empirical success        = ",
        summary.n_closed_loop_success,
        " (",
        round(100 * summary.closed_loop_success_rate; digits = 1),
        "%)",
    )
    println(
        "final in target          = ",
        summary.n_final_in_target,
        " (",
        round(100 * summary.final_target_rate; digits = 1),
        "%)",
    )
    println(
        "obstacle hit             = ",
        summary.n_hit_obstacle,
        " (",
        round(100 * summary.obstacle_rate; digits = 1),
        "%)",
    )
    println(
        "domain exit              = ",
        summary.n_left_domain,
        " (",
        round(100 * summary.domain_exit_rate; digits = 1),
        "%)",
    )
    println(
        "input violation          = ",
        summary.n_input_violation,
        " (",
        round(100 * summary.input_violation_rate; digits = 1),
        "%)",
    )
    println(
        "certificate exit         = ",
        summary.n_ellipsoid_exit,
        " (",
        round(100 * summary.ellipsoid_exit_rate; digits = 1),
        "%)",
    )
    println(
        "outside initial ellipsoid = ",
        summary.n_started_outside_initial_ellipsoid,
        " (",
        round(100 * (1.0 - summary.initial_ellipsoid_coverage_rate); digits = 1),
        "%)",
    )
    println("exit before step         = ", summary.n_exit_before_step)
    println("exit after step          = ", summary.n_exit_after_step)
    println(
        "initial ellipsoid coverage = ",
        round(100 * summary.initial_ellipsoid_coverage_rate; digits = 1),
        "%)",
    )
    if hasproperty(stat_result, :sampling_validation)
        sv = stat_result.sampling_validation
        println("initial sampling check:")
        println("  samples                  = ", sv.n_samples)
        if sv.max_ellipsoid_violation !== nothing
            println("  max ellipsoid violation  = ", sv.max_ellipsoid_violation)
            println("  outside ellipsoid        = ", sv.n_outside_ellipsoid)
        end
        if sv.n_outside_initial_set !== nothing
            println("  outside initial set      = ", sv.n_outside_initial_set)
        end
    end

    if !isempty(stat_result.exit_histogram)
        println("first certificate-exit histogram:")
        for k in sort(collect(keys(stat_result.exit_histogram)))
            println("  k = ", k, " -> ", stat_result.exit_histogram[k])
        end
    end

    return nothing
end

"""
    plot_kappa_statistical_rollouts(stat_result; kwargs...)

Plot the sampled rollouts in a given 2D projection. Success rollouts are green,
obstacle hits are red, ellipsoid exits are orange, and target misses are blue.
"""
function plot_kappa_statistical_rollouts(
    stat_result;
    output_dir::AbstractString = pwd(),
    dims = (1, 2),
    filename = nothing,
    title = nothing,
    show_domain::Bool = true,
    show_sets::Bool = true,
    show_ellipsoids::Bool = false,
    show_exit_points::Bool = false,
    max_ellipsoids::Int = 40,
    wrap_angles::Bool = false,
    unwrap_angles::Bool = false,
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
    rollout_lw::Real = 0.45,
    rollout_alpha::Real = 0.18,
    rollout_ms::Real = 1.2,
    show_initial_samples::Bool = true,
    initial_sample_ms::Real = 1.5,
    initial_sample_alpha::Real = 0.48,
    legend_position = :topleft,
    axis_labels = nothing,
    show_title::Bool = false,
    show_counts_in_legend::Bool = true,
)
    isempty(stat_result.x_rollouts) && error("No rollout available to plot.")
    if output_dir == pwd() &&
       hasproperty(stat_result, :outputs) &&
       stat_result.outputs !== nothing &&
       hasproperty(stat_result.outputs, :plots_dir)
        output_dir = stat_result.outputs.plots_dir
    end
    mkpath(output_dir)

    periodic_cfg = _resolve_periodic_plot_config(
        stat_result;
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )

    if wrap_angles && periodic_cfg.periodic_dims === nothing
        error("wrap_angles = true requires periodic data.")
    end

    d1, d2 = dims
    raw_reference_states =
        hasproperty(stat_result, :reference_states) &&
        !isempty(stat_result.reference_states) ? stat_result.reference_states : nothing
    reference_plot_states =
        (raw_reference_states !== nothing && (unwrap_angles || wrap_angles)) ?
        _unwrap_states_for_angle_plot(
            raw_reference_states,
            periodic_cfg.periodic_dims,
            periodic_cfg.periodic_periods,
        ) : nothing
    should_lift_to_reference =
        reference_plot_states !== nothing &&
        _projection_has_periodic_dim(
            dims,
            periodic_cfg.periodic_dims,
            periodic_cfg.periodic_periods,
        ) &&
        (unwrap_angles || wrap_angles)
    title_txt =
        title !== nothing ? title : (show_title ? "Empirical closed-loop validation" : "")
    xlabel_txt, ylabel_txt =
        axis_labels === nothing ? (_state_axis_label(d1), _state_axis_label(d2)) :
        axis_labels

    fig = plot(;
        thesis_plot_kwargs(; legend = legend_position, size = (760, 560))...,
        aspect_ratio = :equal,
        title = title_txt,
        xlabel = xlabel_txt,
        ylabel = ylabel_txt,
    )

    raw_domain =
        show_domain && hasproperty(stat_result.problem.system, :X) ?
        stat_result.problem.system.X : nothing
    raw_initial_set =
        show_sets && hasproperty(stat_result.problem, :initial_set) ?
        stat_result.problem.initial_set : nothing
    raw_target_set =
        show_sets && hasproperty(stat_result.problem, :target_set) ?
        stat_result.problem.target_set : nothing
    raw_obstacle_set =
        show_sets && hasproperty(stat_result, :obstacle_set) ? stat_result.obstacle_set :
        nothing

    domain = maybe_wrap_set(
        raw_domain,
        periodic_cfg.periodic_dims,
        periodic_cfg.periodic_periods,
        periodic_cfg.periodic_start;
        enabled = wrap_angles,
    )
    initial_set = maybe_wrap_set(
        raw_initial_set,
        periodic_cfg.periodic_dims,
        periodic_cfg.periodic_periods,
        periodic_cfg.periodic_start;
        enabled = wrap_angles,
    )
    initial_set = _lift_set_for_projection(
        initial_set,
        dims;
        ref_state = _reference_state_for_plot(reference_plot_states, 1),
        periodic_dims = periodic_cfg.periodic_dims,
        periodic_periods = periodic_cfg.periodic_periods,
    )
    target_set = maybe_wrap_set(
        raw_target_set,
        periodic_cfg.periodic_dims,
        periodic_cfg.periodic_periods,
        periodic_cfg.periodic_start;
        enabled = wrap_angles,
    )
    target_set = _lift_set_for_projection(
        target_set,
        dims;
        ref_state = _reference_state_for_plot(
            reference_plot_states,
            reference_plot_states === nothing ? 1 : length(reference_plot_states),
        ),
        periodic_dims = periodic_cfg.periodic_dims,
        periodic_periods = periodic_cfg.periodic_periods,
    )
    obstacle_set = maybe_wrap_set(
        raw_obstacle_set,
        periodic_cfg.periodic_dims,
        periodic_cfg.periodic_periods,
        periodic_cfg.periodic_start;
        enabled = wrap_angles,
    )

    _plot_domain_with_obstacles!(fig, domain, [d1, d2])
    initial_set !== nothing && plot!(
        fig,
        initial_set;
        dims = [d1, d2],
        color = PLOT_COLORS[:initial],
        linecolor = PLOT_COLORS[:initial_edge],
        opacity = 0.45,
        lw = 1.0,
        label = L"\mathcal{X}_I",
    )
    target_set !== nothing && plot!(
        fig,
        target_set;
        dims = [d1, d2],
        color = PLOT_COLORS[:target],
        linecolor = PLOT_COLORS[:target_edge],
        opacity = 0.38,
        lw = 1.0,
        label = L"\mathcal{X}_T",
    )
    if obstacle_set !== nothing && !(domain isa UT.LazySetMinus)
        try
            plot!(
                fig,
                obstacle_set;
                dims = [d1, d2],
                color = PLOT_COLORS[:obstacle],
                linecolor = PLOT_COLORS[:obstacle_edge],
                opacity = 0.86,
                lw = 0.9,
                label = L"\mathcal{X}_{obs}",
            )
        catch err
            @debug "Obstacle set could not be plotted." exception = err
        end
    end

    if show_ellipsoids
        ellipsoid_records =
            _ellipsoid_records_for_plot(stat_result; max_keep = max_ellipsoids)

        label_used = false
        for rec in ellipsoid_records
            E =
                should_lift_to_reference ?
                _lift_ellipsoid_for_projection(
                    rec.ellipsoid,
                    dims;
                    ref_state = _reference_state_for_plot(
                        reference_plot_states,
                        rec.state_index,
                    ),
                    periodic_dims = periodic_cfg.periodic_dims,
                    periodic_periods = periodic_cfg.periodic_periods,
                ) : rec.ellipsoid
            E2 = project_ellipsoid_2d(E; dims = dims)
            plot!(
                fig,
                E2;
                color = PLOT_COLORS[:ellipsoid],
                linecolor = PLOT_COLORS[:ellipsoid_edge],
                fillcolor = PLOT_COLORS[:ellipsoid],
                opacity = 0.20,
                lw = 0.8,
                label = label_used ? "" : L"\mathrm{certified\ ellipsoids}",
            )
            label_used = true
        end
    end

    label_counts = Dict{String, Int}()
    for stats in stat_result.rollout_stats
        lbl = _rollout_label(stats)
        label_counts[lbl] = get(label_counts, lbl, 0) + 1
    end

    used_labels = Set{String}()
    for i in eachindex(stat_result.x_rollouts)
        stats = stat_result.rollout_stats[i]
        xs =
            should_lift_to_reference ?
            _lift_state_list_to_plot_reference(
                stat_result.x_rollouts[i],
                reference_plot_states,
                periodic_cfg.periodic_dims,
                periodic_cfg.periodic_periods,
            ) :
            _transform_state_list_for_plot(
                stat_result.x_rollouts[i];
                unwrap_angles = unwrap_angles,
                wrap_angles = wrap_angles,
                periodic_dims = periodic_cfg.periodic_dims,
                periodic_periods = periodic_cfg.periodic_periods,
                periodic_start = periodic_cfg.periodic_start,
            )

        lbl = _rollout_label(stats)
        base_label = _rollout_pretty_label(lbl)
        label_txt =
            lbl in used_labels || !_show_rollout_label_in_legend(lbl) ? "" :
            (show_counts_in_legend ? "$(base_label) (n=$(label_counts[lbl]))" : base_label)
        plot!(
            fig,
            [x[d1] for x in xs],
            [x[d2] for x in xs];
            ms = rollout_ms,
            markershape = :circle,
            markerstrokewidth = 0,
            arrows = false,
            color = _rollout_color(stats),
            alpha = _rollout_line_alpha(stats, rollout_alpha),
            linealpha = _rollout_line_alpha(stats, rollout_alpha),
            lw = _rollout_linewidth(stats, rollout_lw),
            label = label_txt,
        )
        push!(used_labels, lbl)
    end

    if raw_reference_states !== nothing
        xs_ref =
            reference_plot_states === nothing ?
            _transform_state_list_for_plot(
                raw_reference_states;
                unwrap_angles = unwrap_angles,
                wrap_angles = wrap_angles,
                periodic_dims = periodic_cfg.periodic_dims,
                periodic_periods = periodic_cfg.periodic_periods,
                periodic_start = periodic_cfg.periodic_start,
            ) : reference_plot_states
        plot!(
            fig,
            [x[d1] for x in xs_ref],
            [x[d2] for x in xs_ref];
            color = PLOT_COLORS[:nominal],
            lw = 2.2,
            markershape = :circle,
            ms = 1.8,
            markerstrokewidth = 0,
            alpha = 0.95,
            label = "nominal trajectory",
        )
    end

    xs_success_end = Float64[]
    ys_success_end = Float64[]
    xs_failure_end = Float64[]
    ys_failure_end = Float64[]
    for i in eachindex(stat_result.x_rollouts)
        stats = stat_result.rollout_stats[i]
        xend = if should_lift_to_reference
            ref = _reference_state_for_plot(
                reference_plot_states,
                length(stat_result.x_rollouts[i]),
            )
            SVector{length(last(stat_result.x_rollouts[i])), Float64}(
                _lift_state_to_plot_reference(
                    last(stat_result.x_rollouts[i]),
                    ref,
                    periodic_cfg.periodic_dims,
                    periodic_cfg.periodic_periods,
                ),
            )
        else
            _transform_single_state_for_plot(
                last(stat_result.x_rollouts[i]);
                unwrap_angles = unwrap_angles,
                wrap_angles = wrap_angles,
                periodic_dims = periodic_cfg.periodic_dims,
                periodic_periods = periodic_cfg.periodic_periods,
                periodic_start = periodic_cfg.periodic_start,
            )
        end
        if _rollout_is_success(stats)
            push!(xs_success_end, xend[d1])
            push!(ys_success_end, xend[d2])
        else
            push!(xs_failure_end, xend[d1])
            push!(ys_failure_end, xend[d2])
        end
    end
    !isempty(xs_success_end) && scatter!(
        fig,
        xs_success_end,
        ys_success_end;
        color = PLOT_COLORS[:success],
        markershape = :circle,
        ms = 1.5,
        markerstrokewidth = 0,
        alpha = 0.38,
        label = "",
    )
    !isempty(xs_failure_end) && scatter!(
        fig,
        xs_failure_end,
        ys_failure_end;
        color = PLOT_COLORS[:failure],
        markershape = :diamond,
        ms = 2.2,
        markerstrokewidth = 0,
        alpha = 0.90,
        label = "failed endpoints",
    )

    if show_initial_samples && !isempty(stat_result.x0_samples)
        xs0 = Float64[]
        ys0 = Float64[]
        for x0 in stat_result.x0_samples
            xp = if should_lift_to_reference
                ref = _reference_state_for_plot(reference_plot_states, 1)
                SVector{length(x0), Float64}(
                    _lift_state_to_plot_reference(
                        x0,
                        ref,
                        periodic_cfg.periodic_dims,
                        periodic_cfg.periodic_periods,
                    ),
                )
            else
                _transform_single_state_for_plot(
                    x0;
                    unwrap_angles = unwrap_angles,
                    wrap_angles = wrap_angles,
                    periodic_dims = periodic_cfg.periodic_dims,
                    periodic_periods = periodic_cfg.periodic_periods,
                    periodic_start = periodic_cfg.periodic_start,
                )
            end
            push!(xs0, xp[d1])
            push!(ys0, xp[d2])
        end
        scatter!(
            fig,
            xs0,
            ys0;
            color = PLOT_COLORS[:nominal],
            markershape = :circle,
            ms = initial_sample_ms,
            markerstrokewidth = 0,
            alpha = min(Float64(initial_sample_alpha), 0.42),
            label = L"x_0\ \mathrm{samples}",
        )
    end

    if show_exit_points && !isempty(stat_result.ellipsoid_exit_points)
        xs_exit = Float64[]
        ys_exit = Float64[]
        for rec in stat_result.ellipsoid_exit_points
            xp = if should_lift_to_reference
                ref = _reference_state_for_plot(reference_plot_states, rec.k)
                SVector{length(rec.state), Float64}(
                    _lift_state_to_plot_reference(
                        rec.state,
                        ref,
                        periodic_cfg.periodic_dims,
                        periodic_cfg.periodic_periods,
                    ),
                )
            else
                _transform_single_state_for_plot(
                    rec.state;
                    unwrap_angles = unwrap_angles,
                    wrap_angles = wrap_angles,
                    periodic_dims = periodic_cfg.periodic_dims,
                    periodic_periods = periodic_cfg.periodic_periods,
                    periodic_start = periodic_cfg.periodic_start,
                )
            end
            push!(xs_exit, xp[d1])
            push!(ys_exit, xp[d2])
        end

        scatter!(
            fig,
            xs_exit,
            ys_exit;
            color = PLOT_COLORS[:failure],
            markershape = :xcross,
            ms = 2.0,
            alpha = 0.90,
            label = "First certificate exit",
        )
    end

    out_name =
        filename === nothing ? "kappa_statistical_rollouts_$(d1)$(d2).pdf" :
        String(filename)
    out_path = joinpath(output_dir, out_name)
    savefig(fig, out_path)
    display(fig)
    return out_path
end

function save_initial_sampling_debug_plot!(
    stat_result;
    output_dir::AbstractString = pwd(),
    filename::AbstractString = "initial_sampling_debug_12.pdf",
    dims = (1, 2),
    axis_labels = nothing,
    wrap_angles::Bool = false,
    unwrap_angles::Bool = false,
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
    tol::Real = 1.0e-8,
)
    if output_dir == pwd() &&
       hasproperty(stat_result, :outputs) &&
       stat_result.outputs !== nothing &&
       hasproperty(stat_result.outputs, :plots_dir)
        output_dir = stat_result.outputs.plots_dir
    end
    mkpath(output_dir)

    periodic_cfg = _resolve_periodic_plot_config(
        stat_result;
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )

    d1, d2 = dims
    raw_reference_states =
        hasproperty(stat_result, :reference_states) &&
        !isempty(stat_result.reference_states) ? stat_result.reference_states : nothing
    reference_plot_states =
        (raw_reference_states !== nothing && (unwrap_angles || wrap_angles)) ?
        _unwrap_states_for_angle_plot(
            raw_reference_states,
            periodic_cfg.periodic_dims,
            periodic_cfg.periodic_periods,
        ) : nothing
    should_lift_to_reference =
        reference_plot_states !== nothing &&
        _projection_has_periodic_dim(
            dims,
            periodic_cfg.periodic_dims,
            periodic_cfg.periodic_periods,
        ) &&
        (unwrap_angles || wrap_angles)
    xlabel_txt, ylabel_txt =
        axis_labels === nothing ? (latexstring("x_{$d1}"), latexstring("x_{$d2}")) :
        axis_labels

    fig = plot(;
        aspect_ratio = :equal,
        legend = :topright,
        title = "",
        xlabel = xlabel_txt,
        ylabel = ylabel_txt,
        grid = true,
        gridalpha = 0.12,
        foreground_color_grid = PLOT_COLORS[:grid],
        framestyle = :box,
        tickfontsize = 10,
        guidefontsize = 12,
        legendfontsize = 9,
        dpi = 300,
        background_color = :white,
        foreground_color_axis = :gray25,
        foreground_color_border = :gray25,
        margin = 3 * Plots.mm,
    )

    initial_set = maybe_wrap_set(
        stat_result.problem.initial_set,
        periodic_cfg.periodic_dims,
        periodic_cfg.periodic_periods,
        periodic_cfg.periodic_start;
        enabled = wrap_angles,
    )
    initial_set = _lift_set_for_projection(
        initial_set,
        dims;
        ref_state = _reference_state_for_plot(reference_plot_states, 1),
        periodic_dims = periodic_cfg.periodic_dims,
        periodic_periods = periodic_cfg.periodic_periods,
    )
    initial_set !== nothing && plot!(
        fig,
        initial_set;
        dims = [d1, d2],
        color = PLOT_COLORS[:initial],
        linecolor = PLOT_COLORS[:initial_edge],
        opacity = 0.45,
        lw = 1.0,
        label = L"\mathcal{X}_I",
    )

    if stat_result.initial_ellipsoid !== nothing
        E0 =
            should_lift_to_reference ?
            _lift_ellipsoid_for_projection(
                stat_result.initial_ellipsoid,
                dims;
                ref_state = _reference_state_for_plot(reference_plot_states, 1),
                periodic_dims = periodic_cfg.periodic_dims,
                periodic_periods = periodic_cfg.periodic_periods,
            ) : stat_result.initial_ellipsoid
        plot!(
            fig,
            project_ellipsoid_2d(E0; dims = dims);
            color = PLOT_COLORS[:ellipsoid],
            linecolor = PLOT_COLORS[:ellipsoid_edge],
            fillcolor = PLOT_COLORS[:ellipsoid],
            opacity = 0.20,
            lw = 0.8,
            label = L"E_0",
        )
    end

    xs_ok = Float64[]
    ys_ok = Float64[]
    xs_bad = Float64[]
    ys_bad = Float64[]
    for x0 in stat_result.x0_samples
        violates = _ellipsoid_rho(stat_result.initial_ellipsoid, x0) > 1.0 + tol
        xp = if should_lift_to_reference
            ref = _reference_state_for_plot(reference_plot_states, 1)
            SVector{length(x0), Float64}(
                _lift_state_to_plot_reference(
                    x0,
                    ref,
                    periodic_cfg.periodic_dims,
                    periodic_cfg.periodic_periods,
                ),
            )
        else
            _transform_single_state_for_plot(
                x0;
                unwrap_angles = unwrap_angles,
                wrap_angles = wrap_angles,
                periodic_dims = periodic_cfg.periodic_dims,
                periodic_periods = periodic_cfg.periodic_periods,
                periodic_start = periodic_cfg.periodic_start,
            )
        end
        if violates
            push!(xs_bad, xp[d1])
            push!(ys_bad, xp[d2])
        else
            push!(xs_ok, xp[d1])
            push!(ys_ok, xp[d2])
        end
    end

    scatter!(
        fig,
        xs_ok,
        ys_ok;
        color = PLOT_COLORS[:nominal],
        markershape = :circle,
        ms = 1.7,
        markerstrokewidth = 0,
        alpha = 0.50,
        label = L"x_0\ \mathrm{samples}",
    )
    if !isempty(xs_bad)
        scatter!(
            fig,
            xs_bad,
            ys_bad;
            color = PLOT_COLORS[:failure],
            markershape = :xcross,
            ms = 3.0,
            alpha = 0.85,
            label = L"x_0 \notin E_0",
        )
    end

    out_path = joinpath(output_dir, filename)
    savefig(fig, out_path)
    display(fig)
    return out_path
end

"""
    save_kappa_statistical_plots!(stat_result; kwargs...)

Convenience wrapper that saves the standard `(1,2)` and `(3,4)` projections
when available.
"""
function save_kappa_statistical_plots!(
    stat_result;
    output_dir::AbstractString = pwd(),
    basename::AbstractString = "kappa_statistical_rollouts",
    wrap_angles::Bool = false,
    unwrap_angles::Bool = false,
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
    show_ellipsoids::Bool = false,
    show_exit_points::Bool = false,
    max_ellipsoids::Int = 40,
    axis_labels_12 = nothing,
    axis_labels_34 = nothing,
    rollout_lw::Real = 0.6,
    rollout_alpha::Real = 0.24,
    rollout_ms::Real = 1.2,
    show_initial_samples::Bool = true,
    initial_sample_ms::Real = 1.5,
    initial_sample_alpha::Real = 0.48,
    legend_position = :topleft,
    show_title::Bool = false,
    show_counts_in_legend::Bool = true,
)
    if output_dir == pwd() &&
       hasproperty(stat_result, :outputs) &&
       stat_result.outputs !== nothing &&
       hasproperty(stat_result.outputs, :plots_dir)
        output_dir = stat_result.outputs.plots_dir
    end

    x0 = first(stat_result.x_rollouts[1])
    nx = length(x0)

    path12 = plot_kappa_statistical_rollouts(
        stat_result;
        output_dir = output_dir,
        dims = (1, 2),
        filename = "$(basename)_12.pdf",
        wrap_angles = wrap_angles,
        unwrap_angles = unwrap_angles,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
        show_ellipsoids = show_ellipsoids,
        show_exit_points = show_exit_points,
        max_ellipsoids = max_ellipsoids,
        axis_labels = axis_labels_12,
        rollout_lw = rollout_lw,
        rollout_alpha = rollout_alpha,
        rollout_ms = rollout_ms,
        show_initial_samples = show_initial_samples,
        initial_sample_ms = initial_sample_ms,
        initial_sample_alpha = initial_sample_alpha,
        legend_position = legend_position,
        show_title = show_title,
        show_counts_in_legend = show_counts_in_legend,
    )

    path34 = nothing
    if nx >= 4
        path34 = plot_kappa_statistical_rollouts(
            stat_result;
            output_dir = output_dir,
            dims = (3, 4),
            filename = "$(basename)_34.pdf",
            wrap_angles = wrap_angles,
            unwrap_angles = unwrap_angles,
            periodic_dims = periodic_dims,
            periodic_periods = periodic_periods,
            periodic_start = periodic_start,
            show_ellipsoids = show_ellipsoids,
            show_exit_points = show_exit_points,
            max_ellipsoids = max_ellipsoids,
            axis_labels = axis_labels_34,
            rollout_lw = rollout_lw,
            rollout_alpha = rollout_alpha,
            rollout_ms = rollout_ms,
            show_initial_samples = show_initial_samples,
            initial_sample_ms = initial_sample_ms,
            initial_sample_alpha = initial_sample_alpha,
            legend_position = legend_position,
            show_title = show_title,
            show_counts_in_legend = show_counts_in_legend,
        )
    end

    return (; path12, path34)
end
