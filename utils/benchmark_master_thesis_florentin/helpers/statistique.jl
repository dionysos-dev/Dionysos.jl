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

function _validate_initial_samples(sample_set, pts; initial_set = nothing, tol::Real = 1.0e-8)
    if sample_set isa UT.Ellipsoid
        return _validate_ellipsoid_samples(sample_set, pts; tol = tol, initial_set = initial_set)
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
    set isa UT.HyperRectangle && return _sample_points_uniform_in_hyperrectangle(set, n; rng = rng)
    set isa UT.Ellipsoid && return _sample_points_uniform_in_ellipsoid(set, n; rng = rng)

    error("Unsupported sampling set type: $(typeof(set)).")
end

function _project_input_to_domain(u::AbstractVector, u_domain)
    u_domain isa UT.HyperRectangle || return collect(Float64, u)
    return collect(clamp.(u, u_domain.lb, u_domain.ub))
end

function _build_certified_kappa_chain(cert_result)
    cert_result === nothing && error("No certification result available.")
    hasproperty(cert_result, :steps) || error("Certification result has no `steps` field.")

    valid_steps = [
        step for step in cert_result.steps if
        step.status == :ok && step.kappa !== nothing && step.ellipsoid !== nothing
    ]
    isempty(valid_steps) && error("No certified kappa controller found in the certification result.")

    sort!(valid_steps; by = step -> step.k)
    k_sequence = [step.k for step in valid_steps]
    expected_k = collect(first(k_sequence):last(k_sequence))
    k_sequence == expected_k ||
        error("The kappa chain is not contiguous: got $(k_sequence), expected $(expected_k).")
    first(k_sequence) == 1 ||
        error("The kappa chain does not start at k = 1, so it cannot be replayed from the initial set.")

    hasproperty(cert_result, :lmi_data) &&
        cert_result.lmi_data !== nothing &&
        hasproperty(cert_result.lmi_data, :ellipsoids) &&
        !isempty(cert_result.lmi_data.ellipsoids) ||
        error("The certification result does not expose its ellipsoid chain in `lmi_data.ellipsoids`.")

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

function _build_periodic_maps(; periodic_dims = nothing, periodic_periods = nothing, periodic_start = nothing)
    if periodic_dims === nothing || periodic_periods === nothing
        wrap_state = x -> collect(Float64, x)
        lift_state_near_reference = (x, _) -> collect(Float64, x)
        return (; has_periodic = false, wrap_state, lift_state_near_reference)
    end

    start = periodic_start === nothing ? zeros(length(periodic_dims)) : collect(Float64, periodic_start)

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
        disc = ST.discretize_continuous_system(system, Float64(Ts); num_substeps = num_substeps)
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
    stats.first_exit_k !== nothing && return "ellipsoid_exit"
    stats.started_outside_initial_ellipsoid && return "starts_outside"
    stats.final_in_target || return "miss_target"
    return "other"
end

function _rollout_color(stats)
    stats.certified_success && return :seagreen
    stats.closed_loop_success && return :steelblue
    stats.hit_obstacle && return :crimson
    stats.left_domain && return :firebrick
    stats.first_exit_k !== nothing && return :darkorange
    stats.started_outside_initial_ellipsoid && return :goldenrod
    stats.final_in_target || return :mediumpurple
    return :grey45
end

function _rollout_pretty_label(lbl::String)
    lbl == "certified_success" && return "Certified success"
    lbl == "target_reached" && return "Empirical success"
    lbl == "obstacle" && return "Obstacle hit"
    lbl == "domain_exit" && return "Domain exit"
    lbl == "ellipsoid_exit" && return "Certificate exit"
    lbl == "starts_outside" && return "Starts outside initial ellipsoid"
    lbl == "miss_target" && return "Target missed"
    return "Other"
end

function _resolve_periodic_plot_config(
    stat_result;
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
)
    cfg =
        hasproperty(stat_result, :config) ? getproperty(stat_result, :config) :
        (hasproperty(stat_result, :run_config) ? getproperty(stat_result, :run_config) : nothing)

    if periodic_dims === nothing && cfg !== nothing && hasproperty(cfg, :periodic_dims)
        periodic_dims = cfg.periodic_dims
    end
    if periodic_periods === nothing && cfg !== nothing && hasproperty(cfg, :periodic_periods)
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
    num_substeps::Int = 1,
    project_inputs::Bool = true,
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
    n_samples >= 1 || error("n_samples must be >= 1.")
    num_substeps >= 1 || error("num_substeps must be >= 1.")

    hasproperty(run_result, :problem) || error("`run_result` must expose a `problem` field.")
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

    maps = _build_periodic_maps(
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    f_disc = _build_discrete_system_map(problem.system, cert_candidate.Ts; num_substeps = num_substeps)
    x0_samples = sample_points_uniform_in_set(sample_set, n_samples; rng = rng)
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

            if x_next ∈ problem.target_set
                target_reached = true
                target_reach_k = k + 1
                x = x_next
                stop_on_target && break
            end

            x = x_next
        end

        final_state = maps.wrap_state(x)
        final_in_target = final_state ∈ problem.target_set
        target_reached |= final_in_target
        closed_loop_success = !hit_obstacle && !left_domain && target_reached
        certified_success =
            !hit_obstacle &&
            !left_domain &&
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
        final_target_rate = n_target / n_samples,
        ellipsoid_exit_rate = n_ellipsoid_exit / n_samples,
        initial_ellipsoid_coverage_rate = 1.0 - n_started_outside_initial_ellipsoid / n_samples,
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
    rollout_lw::Real = 0.6,
    rollout_alpha::Real = 0.24,
    rollout_ms::Real = 1.2,
    show_initial_samples::Bool = true,
    initial_sample_ms::Real = 1.5,
    initial_sample_alpha::Real = 0.48,
    legend_position = :topright,
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
    title_txt =
        title !== nothing ? title :
        (show_title ? "Empirical closed-loop validation" : "")
    xlabel_txt, ylabel_txt =
        axis_labels === nothing ? (latexstring("x_{$d1}"), latexstring("x_{$d2}")) : axis_labels

    fig = plot(
        ;
        aspect_ratio = :equal,
        legend = legend_position,
        title = title_txt,
        xlabel = xlabel_txt,
        ylabel = ylabel_txt,
        grid = true,
        gridalpha = 0.12,
        foreground_color_grid = :gray82,
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

    raw_domain =
        show_domain && hasproperty(stat_result.problem.system, :X) ? stat_result.problem.system.X : nothing
    raw_initial_set =
        show_sets && hasproperty(stat_result.problem, :initial_set) ? stat_result.problem.initial_set : nothing
    raw_target_set =
        show_sets && hasproperty(stat_result.problem, :target_set) ? stat_result.problem.target_set : nothing
    raw_obstacle_set =
        show_sets && hasproperty(stat_result, :obstacle_set) ? stat_result.obstacle_set : nothing

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
    target_set = maybe_wrap_set(
        raw_target_set,
        periodic_cfg.periodic_dims,
        periodic_cfg.periodic_periods,
        periodic_cfg.periodic_start;
        enabled = wrap_angles,
    )
    obstacle_set = maybe_wrap_set(
        raw_obstacle_set,
        periodic_cfg.periodic_dims,
        periodic_cfg.periodic_periods,
        periodic_cfg.periodic_start;
        enabled = wrap_angles,
    )

    domain !== nothing &&
        plot!(fig, domain; dims = [d1, d2], color = :grey, opacity = 0.08, lw = 0.6, label = L"\mathcal{X}")
    initial_set !== nothing &&
        plot!(fig, initial_set; dims = [d1, d2], color = :seagreen, opacity = 0.14, lw = 0.7, label = L"\mathcal{X}_I")
    target_set !== nothing &&
        plot!(fig, target_set; dims = [d1, d2], color = :steelblue, opacity = 0.18, lw = 0.7, label = L"\mathcal{X}_T")
    if obstacle_set !== nothing
        try
            plot!(
                fig,
                obstacle_set;
                dims = [d1, d2],
                color = :crimson,
                opacity = 0.16,
                lw = 0.7,
                label = L"\mathcal{X}_{obs}",
            )
        catch err
            @debug "Obstacle set could not be plotted." exception = err
        end
    end

    if show_ellipsoids
        ells = _ellipsoids_for_plot(stat_result; max_keep = max_ellipsoids)
        if !isempty(ells) && unwrap_angles
            ells = unwrap_ellipsoid_centers(ells, periodic_cfg.periodic_dims, periodic_cfg.periodic_periods)
        end
        if !isempty(ells) && wrap_angles
            ells = wrap_ellipsoid_centers(
                ells,
                periodic_cfg.periodic_dims,
                periodic_cfg.periodic_periods,
                periodic_cfg.periodic_start,
            )
        end

        label_used = false
        for E in ells
            E2 = project_ellipsoid_2d(E; dims = dims)
            plot!(
                fig,
                E2;
                color = :orange,
                opacity = 0.10,
                lw = 0.6,
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
        xs = _transform_state_list_for_plot(
            stat_result.x_rollouts[i];
            unwrap_angles = unwrap_angles,
            wrap_angles = wrap_angles,
            periodic_dims = periodic_cfg.periodic_dims,
            periodic_periods = periodic_cfg.periodic_periods,
            periodic_start = periodic_cfg.periodic_start,
        )

        lbl = _rollout_label(stats)
        base_label = _rollout_pretty_label(lbl)
        label_txt = lbl in used_labels ? "" :
            (show_counts_in_legend ? "$(base_label) (n=$(label_counts[lbl]))" : base_label)
        plot!(
            fig,
            [x[d1] for x in xs],
            [x[d2] for x in xs];
            lw = rollout_lw,
            ms = rollout_ms,
            markershape = :circle,
            markerstrokewidth = 0,
            arrows = false,
            color = _rollout_color(stats),
            alpha = rollout_alpha,
            label = label_txt,
        )
        push!(used_labels, lbl)
    end

    if show_initial_samples && !isempty(stat_result.x0_samples)
        xs0 = Float64[]
        ys0 = Float64[]
        for x0 in stat_result.x0_samples
            xp = _transform_single_state_for_plot(
                x0;
                unwrap_angles = unwrap_angles,
                wrap_angles = wrap_angles,
                periodic_dims = periodic_cfg.periodic_dims,
                periodic_periods = periodic_cfg.periodic_periods,
                periodic_start = periodic_cfg.periodic_start,
            )
            push!(xs0, xp[d1])
            push!(ys0, xp[d2])
        end
        scatter!(
            fig,
            xs0,
            ys0;
            color = :black,
            markershape = :circle,
            ms = initial_sample_ms,
            markerstrokewidth = 0,
            alpha = initial_sample_alpha,
            label = L"x_0\ \mathrm{samples}",
        )
    end

    if show_exit_points && !isempty(stat_result.ellipsoid_exit_points)
        xs_exit = Float64[]
        ys_exit = Float64[]
        for rec in stat_result.ellipsoid_exit_points
            xp = _transform_single_state_for_plot(
                rec.state;
                unwrap_angles = unwrap_angles,
                wrap_angles = wrap_angles,
                periodic_dims = periodic_cfg.periodic_dims,
                periodic_periods = periodic_cfg.periodic_periods,
                periodic_start = periodic_cfg.periodic_start,
            )
            push!(xs_exit, xp[d1])
            push!(ys_exit, xp[d2])
        end

        scatter!(
            fig,
            xs_exit,
            ys_exit;
            color = :black,
            markershape = :xcross,
            ms = 2.0,
            alpha = 0.55,
            label = "First certificate exit",
        )
    end

    out_name =
        filename === nothing ? "kappa_statistical_rollouts_$(d1)$(d2).pdf" : String(filename)
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
    xlabel_txt, ylabel_txt =
        axis_labels === nothing ? (latexstring("x_{$d1}"), latexstring("x_{$d2}")) : axis_labels

    fig = plot(
        ;
        aspect_ratio = :equal,
        legend = :topright,
        title = "",
        xlabel = xlabel_txt,
        ylabel = ylabel_txt,
        grid = true,
        gridalpha = 0.12,
        foreground_color_grid = :gray82,
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
    initial_set !== nothing &&
        plot!(fig, initial_set; dims = [d1, d2], color = :seagreen, opacity = 0.14, lw = 0.7, label = L"\mathcal{X}_I")

    if stat_result.initial_ellipsoid !== nothing
        E0 = stat_result.initial_ellipsoid
        ells = [E0]
        if unwrap_angles
            ells = unwrap_ellipsoid_centers(ells, periodic_cfg.periodic_dims, periodic_cfg.periodic_periods)
        end
        if wrap_angles
            ells = wrap_ellipsoid_centers(
                ells,
                periodic_cfg.periodic_dims,
                periodic_cfg.periodic_periods,
                periodic_cfg.periodic_start,
            )
        end
        plot!(
            fig,
            project_ellipsoid_2d(first(ells); dims = dims);
            color = :darkorange,
            opacity = 0.12,
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
        xp = _transform_single_state_for_plot(
            x0;
            unwrap_angles = unwrap_angles,
            wrap_angles = wrap_angles,
            periodic_dims = periodic_cfg.periodic_dims,
            periodic_periods = periodic_cfg.periodic_periods,
            periodic_start = periodic_cfg.periodic_start,
        )
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
        color = :black,
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
            color = :crimson,
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
    legend_position = :topright,
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
