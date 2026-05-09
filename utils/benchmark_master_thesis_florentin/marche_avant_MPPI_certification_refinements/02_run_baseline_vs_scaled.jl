include(joinpath(@__DIR__, "run_marche_avant.jl"))

using CSV
using DataFrames
using JLD2
using LinearAlgebra
using Random
using SpecialFunctions: erf
using Statistics
using StaticArrays: SVector

const RESULTS_DIR = joinpath(@__DIR__, "results")
const TRAJECTORY_PATH = joinpath(RESULTS_DIR, "saved_mppi_trajectory.jld2")
const SUMMARY_PATH = joinpath(RESULTS_DIR, "baseline_vs_scaled_summary.csv")
const DETAILED_PATH = joinpath(RESULTS_DIR, "baseline_vs_scaled_detailed.csv")
const BASELINE_RESULT_PATH = joinpath(RESULTS_DIR, "baseline_result.jld2")
const SCALED_RESULT_PATH = joinpath(RESULTS_DIR, "scaled_result.jld2")
const STRESS_SUMMARY_PATH = joinpath(RESULTS_DIR, "statistical_stress_test_summary.csv")
const STRESS_SAMPLES_PATH = joinpath(RESULTS_DIR, "statistical_stress_test_samples.csv")
const STAT_COMPARISON_PATH = joinpath(RESULTS_DIR, "statistical_comparison.csv")

# Set this to a concrete vector, for example [5.0, 5.0, 1.0, 1.0], to override.
const scaled_state_scaling_override = nothing
const n_samples = 1000
const rng_seed = 1234

_vec_string(x) = x === nothing ? "nothing" : join(string.(Float64.(x)), ";")

function _trajectory_matrix(traj)
    elems = collect(ST.enum_elems(traj))
    return reduce(hcat, collect.(Float64, elems))
end

function _candidate_from_matrices(x_traj, u_traj, Ts)
    nx = size(x_traj, 1)
    nu = size(u_traj, 1)
    xs = [SVector{nx, Float64}(x_traj[:, k]) for k in axes(x_traj, 2)]
    us = [SVector{nu, Float64}(u_traj[:, k]) for k in axes(u_traj, 2)]
    return OP.CandidateTrajectory(
        ST.Trajectory(xs),
        ST.Trajectory(us);
        Ts = Ts,
        source = :saved_mppi,
        metadata = (; saved_from = TRAJECTORY_PATH),
    )
end

function _load_candidate()
    data = load(TRAJECTORY_PATH)
    cand = _candidate_from_matrices(data["x_traj"], data["u_traj"], data["Ts"])
    return default_prepare_for_certification(MarcheAvantConfig())(cand), data
end

function _default_scaled_state_scaling(x_traj)
    scaled_state_scaling_override !== nothing && return Float64.(scaled_state_scaling_override)
    n = size(x_traj, 1)
    s = ones(n)
    if n >= 1
        s[1] = max(maximum(x_traj[1, :]) - minimum(x_traj[1, :]), 1e-2)
    end
    if n >= 2
        s[2] = max(maximum(x_traj[2, :]) - minimum(x_traj[2, :]), 1e-2)
    end
    return s
end

function _build_certifier(problem, system_cfg, cfg; state_scaling)
    opts = (
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        ΔX = cfg.ΔX,
        ΔU = cfg.ΔU,
        ΔW = cfg.ΔW,
        rayon_terminal = cfg.terminal_radius,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
        state_scaling = state_scaling === nothing ? nothing : Float64.(state_scaling),
    )

    Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    cert_cfg = SC.EllipsoidalBackwardConfig(
        problem.system.X,
        problem.system.U,
        Wdom,
        build_backend(; verbose = cfg.verbose),
        opts,
    )

    return SC.EllipsoidalBackwardCertifier(
        cert_cfg,
        build_symbolic_builder(AV, system_cfg.params),
    )
end

function _assert_close(a, b; atol = 1e-8, rtol = 1e-8, label = "values")
    err = norm(collect(a) .- collect(b))
    scale = max(norm(collect(a)), norm(collect(b)), 1.0)
    err <= atol + rtol * scale || error("Scaling sanity check failed for $(label): error=$(err)")
    return nothing
end

function _run_scaling_sanity_checks(problem, system_cfg, cfg, candidate, state_scaling)
    state_scaling === nothing && return nothing

    cert = _build_certifier(problem, system_cfg, cfg; state_scaling = state_scaling)
    ctx = SC.build_symbolic_context(problem, candidate, cert.config, cert.symbolic_builder)

    xk = collect(ctx.xs[1])
    xnext = collect(ctx.xs[2])
    uk = collect(ctx.us[1])
    wk = zeros(length(ctx.symbolic.w))
    Xbar = IA.IntervalBox(xk .+ ctx.symbolic.ΔX)
    Ubar = IA.IntervalBox(uk .+ ctx.symbolic.ΔU)
    Wbar = IA.IntervalBox(wk .+ ctx.symbolic.ΔW)

    affsys, _ = ST.buildAffineApproximation(
        ctx.symbolic.fsymbolic,
        ctx.symbolic.x,
        ctx.symbolic.u,
        ctx.symbolic.w,
        xk,
        uk,
        wk,
        Xbar,
        Ubar,
        Wbar,
    )

    D = Diagonal(state_scaling)
    Sinv = Diagonal(1.0 ./ state_scaling)
    Az = Sinv * affsys.A * D
    Bz = Sinv * affsys.B
    cz = Sinv * (affsys.A * xk + affsys.c - xnext)
    Dwz = Sinv * affsys.D

    rng = Random.MersenneTwister(42)
    z = randn(rng, length(xk))
    u = uk .+ 0.1 .* randn(rng, length(uk))
    w = zeros(size(affsys.D, 2))
    x = xk + D * z
    lhs = Sinv * (affsys.A * x + affsys.B * u + affsys.c + affsys.D * w - xnext)
    rhs = Az * z + Bz * u + cz + Dwz * w
    _assert_close(lhs, rhs; label = "scaled affine dynamics")

    Pz = Matrix{Float64}(I, length(xk), length(xk))
    Ez = UT.Ellipsoid(Pz, zeros(length(xk)))
    Ex = UT.Ellipsoid(Matrix(Sinv' * Pz * Sinv), xk)
    z_sample = 0.25 .* randn(rng, length(xk))
    x_sample = xk + D * z_sample
    _assert_close(
        [z_sample' * Pz * z_sample],
        [(x_sample - xk)' * Matrix(Ex.P) * (x_sample - xk)];
        label = "ellipsoid round trip",
    )

    nu = length(uk)
    Kz = randn(rng, nu, length(xk))
    ell = randn(rng, nu)
    Kx = Kz * Sinv
    b = ell - Kx * xk
    _assert_close(ell + Kz * z_sample, Kx * x_sample + b; label = "feedback round trip")

    println("scaling sanity checks passed")
    return nothing
end

function _run_certification(method, problem, system_cfg, cfg, candidate; state_scaling)
    cert = _build_certifier(problem, system_cfg, cfg; state_scaling = state_scaling)
    SC.set_problem!(cert, problem)
    SC.set_trajectory!(cert, candidate)
    t0 = time()
    SC.certify!(cert)
    total_time = time() - t0
    result = SC.get_result(cert)
    return (; method, state_scaling, cert, result, total_time)
end

function _volumes_by_k(result)
    rows = NamedTuple[]
    for step in sort(result.steps; by = s -> s.k)
        volume = step.ellipsoid === nothing ? missing : Float64(UT.get_volume(step.ellipsoid))
        solve_time =
            hasproperty(step.summary, :solve_time) && step.summary.solve_time !== nothing ?
            Float64(step.summary.solve_time) : missing
        status =
            hasproperty(step.summary, :solver_status) ? string(step.summary.solver_status) :
            string(step.status)
        push!(rows, (; k = step.k, volume, solve_time, status))
    end
    return rows
end

function _solver_status_summary(result)
    statuses = [
        hasproperty(step.summary, :solver_status) ? string(step.summary.solver_status) :
        string(step.status) for step in result.steps
    ]
    return join(["$(s):$(count(==(s), statuses))" for s in sort(unique(statuses))], ";")
end

function _summary_row(run; horizon)
    rows = _volumes_by_k(run.result)
    vols = [Float64(row.volume) for row in rows if !ismissing(row.volume)]
    solve_times = [Float64(row.solve_time) for row in rows if !ismissing(row.solve_time)]
    return (;
        method = run.method,
        state_scaling = _vec_string(run.state_scaling),
        success = run.result.success,
        number_of_certified_transitions = count(step -> step.status == :ok, run.result.steps),
        horizon,
        min_volume = isempty(vols) ? missing : minimum(vols),
        median_volume = isempty(vols) ? missing : median(vols),
        mean_volume = isempty(vols) ? missing : mean(vols),
        max_volume = isempty(vols) ? missing : maximum(vols),
        initial_volume = isempty(vols) ? missing : first(vols),
        final_volume =
            run.result.lmi_data !== nothing && hasproperty(run.result.lmi_data, :ellipsoids) ?
            Float64(UT.get_volume(first(run.result.lmi_data.ellipsoids))) : missing,
        total_certification_time = run.total_time,
        average_step_time = isempty(solve_times) ? missing : mean(solve_times),
        solver_status_summary = _solver_status_summary(run.result),
    )
end

function _detailed_rows(run)
    return [
        (;
            method = run.method,
            k = row.k,
            volume_Ek = row.volume,
            solver_status = row.status,
            solve_time = row.solve_time,
        ) for row in _volumes_by_k(run.result)
    ]
end

function _run_result_payload(run, problem, cfg, candidate, outputs)
    return (;
        problem,
        result = (; certification = run.result),
        certification_candidate = candidate,
        config = cfg,
        outputs,
    )
end

function _ellipsoid_value(E, x)
    dx = collect(Float64, x) .- collect(Float64, E.c)
    return Float64(dx' * Matrix(E.P) * dx)
end

function _final_distance_to_target(problem, x)
    hasproperty(problem.target_set, :lb) && hasproperty(problem.target_set, :ub) || return missing
    center = (collect(Float64, problem.target_set.lb) .+ collect(Float64, problem.target_set.ub)) ./ 2
    return norm(collect(Float64, x) .- center)
end

function _simulate_samples(run_payload, x0_samples; method, scenario)
    problem = run_payload.problem
    cert_result = run_payload.result.certification
    chain = _build_certified_kappa_chain(cert_result)
    cert_candidate = run_payload.certification_candidate
    cert_xs = collect(ST.enum_elems(cert_candidate.x_traj))
    cfg = run_payload.config

    maps = _build_periodic_maps(
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
    )
    f_disc = _build_discrete_system_map(
        problem.system,
        cert_candidate.Ts;
        num_substeps = cfg.symbolic_rk4_substeps,
    )
    domain_set = _domain_set(problem.system)
    obstacle_set = _obstacle_set(problem.system)

    rows = Vector{Dict{Symbol, Any}}()
    rollouts = Vector{Vector{Vector{Float64}}}()

    for (sample_id, x0) in enumerate(x0_samples)
        x = maps.wrap_state(x0)
        x_hist = Vector{Vector{Float64}}([copy(x)])
        status = "success"
        first_failure_step = missing
        max_violation = -Inf

        try
            x_init_chain = maps.lift_state_near_reference(x, cert_xs[first(chain.k_sequence)])
            if !(x_init_chain ∈ chain.initial_ellipsoid)
                status = "left_chain"
                first_failure_step = first(chain.k_sequence)
            end

            for k in chain.k_sequence
                x_phys = maps.wrap_state(x)
                x_chain = maps.lift_state_near_reference(x_phys, cert_xs[k])
                max_violation = max(max_violation, _ellipsoid_value(chain.ellipsoid_by_k[k], x_chain))

                if status == "success" && !(x_chain ∈ chain.ellipsoid_by_k[k])
                    status = "left_chain"
                    first_failure_step = k
                end

                u = _eval_kappa_controller(chain.kappa_by_k[k], x_chain)
                hasproperty(problem.system, :U) && (u = _project_input_to_domain(u, problem.system.U))

                x_next = maps.wrap_state(f_disc(x_phys, u))
                push!(x_hist, copy(x_next))

                x_next_chain = maps.lift_state_near_reference(x_next, cert_xs[k + 1])
                max_violation = max(
                    max_violation,
                    _ellipsoid_value(chain.ellipsoid_by_k[k + 1], x_next_chain),
                )

                if status == "success" && !(x_next_chain ∈ chain.ellipsoid_by_k[k + 1])
                    status = "left_chain"
                    first_failure_step = k + 1
                end

                if status == "success" && _point_in_any_set(x_next, obstacle_set)
                    status = "obstacle_violation"
                    first_failure_step = k + 1
                end
                if status == "success" && domain_set !== nothing && !(x_next ∈ domain_set)
                    status = "obstacle_violation"
                    first_failure_step = k + 1
                end

                x = x_next
            end

            final_state = maps.wrap_state(x)
            if status == "success" && !(final_state ∈ problem.target_set)
                status = "target_failed"
                first_failure_step = last(chain.k_sequence) + 1
            end
        catch err
            status = "simulation_error"
            first_failure_step = missing
        end

        final_state = maps.wrap_state(last(x_hist))
        row = Dict{Symbol, Any}(
            :scenario => scenario,
            :method => method,
            :sample_id => sample_id,
            :status => status,
            :first_failure_step => first_failure_step,
            :final_state => join(string.(final_state), ";"),
            :final_distance_to_target => _final_distance_to_target(problem, final_state),
            :maximum_normalized_ellipsoid_violation => max_violation,
        )
        for i in eachindex(x0)
            row[Symbol("x0_", i)] = Float64(x0[i])
        end
        push!(rows, row)
        push!(rollouts, x_hist)
    end

    return (; rows, rollouts)
end

function _wilson_interval(successes, n; z = 1.959963984540054)
    n == 0 && return (missing, missing)
    phat = successes / n
    denom = 1 + z^2 / n
    center = (phat + z^2 / (2n)) / denom
    half = z * sqrt((phat * (1 - phat) + z^2 / (4n)) / n) / denom
    return (max(0.0, center - half), min(1.0, center + half))
end

function _stress_summary(method, sample_rows; scenario)
    n = length(sample_rows)
    count_status(s) = count(row -> row[:status] == s, sample_rows)
    n_success = count_status("success")
    n_left_chain = count_status("left_chain")
    n_target_failed = count_status("target_failed")
    n_obstacle_violation = count_status("obstacle_violation")
    n_simulation_error = count_status("simulation_error")
    success_ci = _wilson_interval(n_success, n)
    left_ci = _wilson_interval(n_left_chain, n)
    return (;
        scenario,
        method,
        n_samples = n,
        n_success,
        n_left_chain,
        n_target_failed,
        n_obstacle_violation,
        n_simulation_error,
        success_rate = n_success / n,
        left_chain_rate = n_left_chain / n,
        success_ci_low = success_ci[1],
        success_ci_high = success_ci[2],
        left_chain_ci_low = left_ci[1],
        left_chain_ci_high = left_ci[2],
    )
end

function _normal_cdf(x)
    return 0.5 * (1 + erf(x / sqrt(2)))
end

function _two_proportion_pvalue(x1, n1, x2, n2)
    (n1 == 0 || n2 == 0) && return missing
    p_pool = (x1 + x2) / (n1 + n2)
    se = sqrt(p_pool * (1 - p_pool) * (1 / n1 + 1 / n2))
    se == 0 && return 1.0
    z = (x2 / n2 - x1 / n1) / se
    return 2 * (1 - _normal_cdf(abs(z)))
end

function _mcnemar_pvalue(baseline_rows, scaled_rows)
    by_id_scaled = Dict(row[:sample_id] => row for row in scaled_rows)
    b_success_s_failure = 0
    b_failure_s_success = 0
    for brow in baseline_rows
        srow = by_id_scaled[brow[:sample_id]]
        bsucc = brow[:status] == "success"
        ssucc = srow[:status] == "success"
        bsucc && !ssucc && (b_success_s_failure += 1)
        !bsucc && ssucc && (b_failure_s_success += 1)
    end
    ndisc = b_success_s_failure + b_failure_s_success
    ndisc == 0 && return 1.0
    stat = (abs(b_success_s_failure - b_failure_s_success) - 1)^2 / ndisc
    return 1 - erf(sqrt(stat / 2))
end

function _comparison_rows(baseline_rows, scaled_rows)
    bsum = _stress_summary("baseline", baseline_rows; scenario = "paired_baseline_E0")
    ssum = _stress_summary("scaled", scaled_rows; scenario = "paired_baseline_E0")
    return [
        (;
            comparison_type = "paired_baseline_E0_mcnemar",
            baseline_success_rate = bsum.success_rate,
            scaled_success_rate = ssum.success_rate,
            difference = ssum.success_rate - bsum.success_rate,
            p_value = _mcnemar_pvalue(baseline_rows, scaled_rows),
            ci_baseline_low = bsum.success_ci_low,
            ci_baseline_high = bsum.success_ci_high,
            ci_scaled_low = ssum.success_ci_low,
            ci_scaled_high = ssum.success_ci_high,
            interpretation = "Paired stress test on samples drawn from baseline E0.",
        ),
        (;
            comparison_type = "two_proportion_success",
            baseline_success_rate = bsum.success_rate,
            scaled_success_rate = ssum.success_rate,
            difference = ssum.success_rate - bsum.success_rate,
            p_value = _two_proportion_pvalue(
                bsum.n_success,
                bsum.n_samples,
                ssum.n_success,
                ssum.n_samples,
            ),
            ci_baseline_low = bsum.success_ci_low,
            ci_baseline_high = bsum.success_ci_high,
            ci_scaled_low = ssum.success_ci_low,
            ci_scaled_high = ssum.success_ci_high,
            interpretation = "Approximate unpaired two-proportion z-test; paired McNemar row is preferred for common samples.",
        ),
    ]
end

function main()
    isfile(TRAJECTORY_PATH) ||
        error("Missing saved trajectory. Run 01_generate_and_save_mppi_trajectory.jl first.")
    mkpath(RESULTS_DIR)

    cfg = MarcheAvantConfig(; output_root = RESULTS_DIR, plot_gif = false)
    candidate, data = _load_candidate()
    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    problem = build_problem(system_cfg, control_cfg)
    outputs = (; root = RESULTS_DIR, plots_dir = joinpath(RESULTS_DIR, "plots"))
    mkpath(outputs.plots_dir)

    scaled_state_scaling = _default_scaled_state_scaling(data["x_traj"])
    println("scaled_state_scaling = ", scaled_state_scaling)
    _run_scaling_sanity_checks(problem, system_cfg, cfg, candidate, scaled_state_scaling)

    baseline = _run_certification(
        "baseline",
        problem,
        system_cfg,
        cfg,
        candidate;
        state_scaling = nothing,
    )
    scaled = _run_certification(
        "scaled",
        problem,
        system_cfg,
        cfg,
        candidate;
        state_scaling = scaled_state_scaling,
    )

    summary = DataFrame([
        _summary_row(baseline; horizon = OP.horizon(candidate)),
        _summary_row(scaled; horizon = OP.horizon(candidate)),
    ])
    detailed = DataFrame(vcat(_detailed_rows(baseline), _detailed_rows(scaled)))
    CSV.write(SUMMARY_PATH, summary)
    CSV.write(DETAILED_PATH, detailed)

    baseline_payload = _run_result_payload(baseline, problem, cfg, candidate, outputs)
    scaled_payload = _run_result_payload(scaled, problem, cfg, candidate, outputs)
    jldsave(BASELINE_RESULT_PATH; run = baseline_payload, certification = baseline.result)
    jldsave(SCALED_RESULT_PATH; run = scaled_payload, certification = scaled.result)

    if baseline.result.success && scaled.result.success
        rng = Random.MersenneTwister(rng_seed)
        baseline_chain = _build_certified_kappa_chain(baseline.result)
        scaled_chain = _build_certified_kappa_chain(scaled.result)

        baseline_x0_samples = sample_points_uniform_in_set(
            baseline_chain.initial_ellipsoid,
            n_samples;
            rng = rng,
        )
        scaled_x0_samples = sample_points_uniform_in_set(
            scaled_chain.initial_ellipsoid,
            n_samples;
            rng = rng,
        )

        baseline_own_stress = _simulate_samples(
            baseline_payload,
            baseline_x0_samples;
            method = "baseline",
            scenario = "own_ellipsoid",
        )
        scaled_own_stress = _simulate_samples(
            scaled_payload,
            scaled_x0_samples;
            method = "scaled",
            scenario = "own_ellipsoid",
        )
        baseline_paired_stress = _simulate_samples(
            baseline_payload,
            baseline_x0_samples;
            method = "baseline",
            scenario = "paired_baseline_E0",
        )
        scaled_paired_stress = _simulate_samples(
            scaled_payload,
            baseline_x0_samples;
            method = "scaled",
            scenario = "paired_baseline_E0",
        )

        stress_rows = vcat(
            baseline_own_stress.rows,
            scaled_own_stress.rows,
            baseline_paired_stress.rows,
            scaled_paired_stress.rows,
        )
        CSV.write(STRESS_SAMPLES_PATH, DataFrame(stress_rows))
        CSV.write(
            STRESS_SUMMARY_PATH,
            DataFrame([
                _stress_summary(
                    "baseline",
                    baseline_own_stress.rows;
                    scenario = "own_ellipsoid",
                ),
                _stress_summary("scaled", scaled_own_stress.rows; scenario = "own_ellipsoid"),
                _stress_summary(
                    "baseline",
                    baseline_paired_stress.rows;
                    scenario = "paired_baseline_E0",
                ),
                _stress_summary(
                    "scaled",
                    scaled_paired_stress.rows;
                    scenario = "paired_baseline_E0",
                ),
            ]),
        )
        CSV.write(
            STAT_COMPARISON_PATH,
            DataFrame(_comparison_rows(baseline_paired_stress.rows, scaled_paired_stress.rows)),
        )
        jldsave(
            joinpath(RESULTS_DIR, "statistical_rollouts.jld2");
            baseline_rollouts = baseline_paired_stress.rollouts,
            scaled_rollouts = scaled_paired_stress.rollouts,
        )
    else
        println("Skipping statistical stress test because one certification failed.")
    end

    println("summary: ", SUMMARY_PATH)
    println("details: ", DETAILED_PATH)
    println("baseline result: ", BASELINE_RESULT_PATH)
    println("scaled result: ", SCALED_RESULT_PATH)
    return (; baseline, scaled)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
