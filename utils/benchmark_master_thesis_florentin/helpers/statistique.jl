import LinearAlgebra as LA
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
)
    n >= 1 || error("n must be >= 1.")

    nx = length(E.c)
    L = LA.cholesky(LA.Symmetric(Matrix{Float64}(E.P))).L
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
        pts[i] = collect(c .+ L \ z)
    end

    return pts
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

    k_terminal = last(k_sequence) + 1
    ellipsoid_by_k[k_terminal] = cert_result.lmi_data.ellipsoids[1]

    return (;
        k_sequence = k_sequence,
        kappa_by_k = kappa_by_k,
        ellipsoid_by_k = ellipsoid_by_k,
        initial_ellipsoid = ellipsoid_by_k[1],
        target_ellipsoid = ellipsoid_by_k[k_terminal],
    )
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
    stats.first_exit_k !== nothing && return "ellipsoid_exit"
    stats.started_outside_initial_ellipsoid && return "starts_outside"
    stats.final_in_target || return "miss_target"
    return "other"
end

function _rollout_color(stats)
    stats.certified_success && return :seagreen3
    stats.closed_loop_success && return :deepskyblue3
    stats.hit_obstacle && return :crimson
    stats.first_exit_k !== nothing && return :darkorange2
    stats.started_outside_initial_ellipsoid && return :goldenrod3
    stats.final_in_target || return :royalblue3
    return :grey50
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
- points are sampled from `problem.initial_set` by default,
- for periodic coordinates, pass `periodic_dims`, `periodic_periods`,
  and `periodic_start` so that the state is wrapped for the concrete
  simulation and lifted near the certification trajectory for ellipsoid
  and kappa evaluation.
"""
function run_kappa_statistical_check(
    run_result;
    n_samples::Int = 100,
    sample_set = nothing,
    sample_from::Symbol = :initial_set,
    num_substeps::Int = 1,
    project_inputs::Bool = true,
    check_domain::Bool = true,
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
        if sample_from == :initial_set
            sample_set = problem.initial_set
        elseif sample_from == :initial_ellipsoid
            sample_set = chain.initial_ellipsoid
        else
            error("Unsupported `sample_from` mode: $(sample_from).")
        end
    end

    maps = _build_periodic_maps(
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
    f_disc = _build_discrete_system_map(problem.system, cert_candidate.Ts; num_substeps = num_substeps)
    x0_samples = sample_points_uniform_in_set(sample_set, n_samples; rng = rng)

    x_rollouts = Vector{Vector{Vector{Float64}}}(undef, n_samples)
    u_rollouts = Vector{Vector{Vector{Float64}}}(undef, n_samples)
    rollout_stats = Vector{NamedTuple}(undef, n_samples)
    ellipsoid_exit_points = NamedTuple[]

    for s in 1:n_samples
        x = maps.wrap_state(x0_samples[s])
        x_hist = Vector{Vector{Float64}}([copy(x)])
        u_hist = Vector{Vector{Float64}}()

        hit_obstacle = false
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

            if check_domain && hasproperty(problem.system, :X) && !(x_next ∈ problem.system.X)
                hit_obstacle = true
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
        closed_loop_success = !hit_obstacle && target_reached
        certified_success =
            !hit_obstacle &&
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
        n_final_in_target = n_target,
        n_ellipsoid_exit = n_ellipsoid_exit,
        n_exit_before_step = n_exit_before,
        n_exit_after_step = n_exit_after,
        n_started_outside_initial_ellipsoid = n_started_outside_initial_ellipsoid,
        n_outside_initial_ellipsoid = n_started_outside_initial_ellipsoid,
        success_rate = n_certified_success / n_samples,
        certified_success_rate = n_certified_success / n_samples,
        closed_loop_success_rate = n_closed_loop_success / n_samples,
        obstacle_rate = n_obstacle / n_samples,
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
        certification_result = cert_result,
        certification_candidate = cert_candidate,
        k_sequence = chain.k_sequence,
        kappa_by_k = chain.kappa_by_k,
        ellipsoid_by_k = chain.ellipsoid_by_k,
        initial_ellipsoid = chain.initial_ellipsoid,
        target_ellipsoid = chain.target_ellipsoid,
        reference_states = cert_xs,
        x0_samples = x0_samples,
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

    println("\n=== Kappa statistical check ===")
    println("samples = ", summary.n_samples)
    println(
        "certified_success = ",
        summary.n_certified_success,
        " (",
        round(100 * summary.certified_success_rate; digits = 1),
        "%)",
    )
    println(
        "closed_loop_success = ",
        summary.n_closed_loop_success,
        " (",
        round(100 * summary.closed_loop_success_rate; digits = 1),
        "%)",
    )
    println(
        "final_in_target = ",
        summary.n_final_in_target,
        " (",
        round(100 * summary.final_target_rate; digits = 1),
        "%)",
    )
    println(
        "hit_obstacle = ",
        summary.n_hit_obstacle,
        " (",
        round(100 * summary.obstacle_rate; digits = 1),
        "%)",
    )
    println(
        "ellipsoid_exit = ",
        summary.n_ellipsoid_exit,
        " (",
        round(100 * summary.ellipsoid_exit_rate; digits = 1),
        "%)",
    )
    println("exit_before_step = ", summary.n_exit_before_step)
    println("exit_after_step = ", summary.n_exit_after_step)
    println("started_outside_initial_ellipsoid = ", summary.n_started_outside_initial_ellipsoid)
    println(
        "initial_ellipsoid_coverage = ",
        round(100 * summary.initial_ellipsoid_coverage_rate; digits = 1),
        "%)",
    )

    if !isempty(stat_result.exit_histogram)
        println("first_exit_histogram:")
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
    show_ellipsoids::Bool = true,
    show_exit_points::Bool = true,
    max_ellipsoids::Int = 40,
    wrap_angles::Bool = false,
    unwrap_angles::Bool = false,
    periodic_dims = nothing,
    periodic_periods = nothing,
    periodic_start = nothing,
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
    succ_rate = stat_result.summary.success_rate
    title_txt =
        title !== nothing ? title :
        "Kappa statistics ($(d1),$(d2)) | success=$(round(100 * succ_rate; digits = 1))%"

    fig = plot(; aspect_ratio = :equal, legend = true, title = title_txt)

    raw_domain =
        show_domain && hasproperty(stat_result.problem.system, :X) ? stat_result.problem.system.X : nothing
    raw_initial_set =
        show_sets && hasproperty(stat_result.problem, :initial_set) ? stat_result.problem.initial_set : nothing
    raw_target_set =
        show_sets && hasproperty(stat_result.problem, :target_set) ? stat_result.problem.target_set : nothing

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

    domain !== nothing && plot!(fig, domain; dims = [d1, d2], color = :grey, opacity = 0.10, label = "Domain")
    initial_set !== nothing &&
        plot!(fig, initial_set; dims = [d1, d2], color = :green, opacity = 0.12, label = "Initial set")
    target_set !== nothing &&
        plot!(fig, target_set; dims = [d1, d2], color = :red, opacity = 0.22, label = "Target set")

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
                opacity = 0.12,
                lw = 0.8,
                label = label_used ? "" : "Certified ellipsoids",
            )
            label_used = true
        end
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

        traj = ST.Trajectory(xs)
        lbl = _rollout_label(stats)
        label_txt =
            lbl in used_labels ? "" :
            (lbl == "certified_success" ? "Certified success" :
             lbl == "target_reached" ? "Target reached" :
             lbl == "obstacle" ? "Obstacle hit" :
             lbl == "ellipsoid_exit" ? "Ellipsoid exit" :
             lbl == "starts_outside" ? "Starts outside E1" :
             lbl == "miss_target" ? "Miss target" : "Other")
        plot!(
            fig,
            traj;
            dims = [d1, d2],
            ms = 1.0,
            arrows = false,
            color = _rollout_color(stats),
            alpha = 0.35,
            label = label_txt,
        )
        push!(used_labels, lbl)
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
            marker = :xcross,
            ms = 2.5,
            alpha = 0.75,
            label = "First ellipsoid exit",
        )
    end

    out_name =
        filename === nothing ? "kappa_statistical_rollouts_$(d1)$(d2).pdf" : String(filename)
    out_path = joinpath(output_dir, out_name)
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
    show_ellipsoids::Bool = true,
    show_exit_points::Bool = true,
    max_ellipsoids::Int = 40,
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
        )
    end

    return (; path12, path34)
end
