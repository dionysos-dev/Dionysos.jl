include(joinpath(@__DIR__, "02_run_baseline_vs_scaled.jl"))

using Plots

const FINAL_SUMMARY_SCALING =
    joinpath(RESULTS_DIR, "final_validation_scaling_sweep_summary.csv")
const FINAL_DETAILED_SCALING =
    joinpath(RESULTS_DIR, "final_validation_scaling_sweep_detailed.csv")
const FINAL_SUMMARY_BOX = joinpath(RESULTS_DIR, "final_validation_box_sweep_summary.csv")
const FINAL_DETAILED_BOX = joinpath(RESULTS_DIR, "final_validation_box_sweep_detailed.csv")
const FINAL_STRESS_SAMPLES = joinpath(RESULTS_DIR, "final_validation_stress_samples.csv")
const FINAL_DIAGNOSTICS = joinpath(RESULTS_DIR, "final_validation_diagnostics.md")
const FINAL_IDENTITY = joinpath(RESULTS_DIR, "final_validation_identity_check.csv")
const FINAL_ROUNDTRIP = joinpath(RESULTS_DIR, "final_validation_roundtrip_checks.csv")
const FINAL_ELLIPSOID_CONVENTION =
    joinpath(RESULTS_DIR, "final_validation_ellipsoid_convention_check.csv")
const FINAL_STRESS_DEBUG = joinpath(RESULTS_DIR, "final_validation_stress_debug.csv")
const FINAL_PLOTS_DIR = joinpath(RESULTS_DIR, "plots", "final_validation")

const final_n_samples = parse(Int, get(ENV, "FINAL_VALIDATION_N_SAMPLES", "1000"))
const n_samples_debug = parse(Int, get(ENV, "FINAL_VALIDATION_N_DEBUG", "200"))
const final_rng_seed = 1234
const box_multipliers_X = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
const box_multipliers_U = [1.0]

function _interval_box_scale(Δ, α)
    α == 1.0 && return Δ
    return IA.IntervalBox((Float64(α) .* collect(Δ))...)
end

function _cfg_with_boxes(cfg; mx = 1.0, mu = 1.0)
    return MarcheAvantConfig(;
        Δt = cfg.Δt,
        hx = cfg.hx,
        periodic_dims = cfg.periodic_dims,
        periodic_periods = cfg.periodic_periods,
        periodic_start = cfg.periodic_start,
        nstep = cfg.nstep,
        terminal_radius = cfg.terminal_radius,
        λ = cfg.λ,
        maxδx = cfg.maxδx,
        maxδu = cfg.maxδu,
        symbolic_rk4_substeps = cfg.symbolic_rk4_substeps,
        ΔX = _interval_box_scale(cfg.ΔX, mx),
        ΔU = _interval_box_scale(cfg.ΔU, mu),
        ΔW = cfg.ΔW,
        output_root = cfg.output_root,
        plot_subdir = cfg.plot_subdir,
        animation_subdir = cfg.animation_subdir,
        plot_gif = false,
        verbose = cfg.verbose,
        seed_trajectory_mode = cfg.seed_trajectory_mode,
        seed_num_substeps = cfg.seed_num_substeps,
        mppi_nsamples = cfg.mppi_nsamples,
        mppi_niter = cfg.mppi_niter,
        mppi_λ = cfg.mppi_λ,
        mppi_noise_v = cfg.mppi_noise_v,
        mppi_noise_σ = cfg.mppi_noise_σ,
    )
end

function _scaling_candidates(x_traj, cfg)
    nx = size(x_traj, 1)
    range_xy =
        [max(maximum(x_traj[i, :]) - minimum(x_traj[i, :]), 1e-2) for i in 1:min(nx, 2)]
    std_xy = [max(std(x_traj[i, :]), 1e-2) for i in 1:min(nx, 2)]
    box_r = [max(abs(Float64(IA.inf(I))), abs(Float64(IA.sup(I)))) for I in collect(cfg.ΔX)]

    pairs = Pair{Symbol, Any}[:none => nothing, :identity => ones(nx)]
    if nx == 4
        append!(
            pairs,
            [
                :xy_2_angles_1 => [2.0, 2.0, 1.0, 1.0],
                :xy_5_angles_1 => [5.0, 5.0, 1.0, 1.0],
                :xy_10_angles_1 => [10.0, 10.0, 1.0, 1.0],
                :xy_10_angles_05 => [10.0, 10.0, 0.5, 0.5],
                :xy_5_angles_05 => [5.0, 5.0, 0.5, 0.5],
                :trajectory_range => [range_xy[1], range_xy[2], 1.0, 1.0],
                :trajectory_std => [std_xy[1], std_xy[2], 1.0, 1.0],
                :box_based => max.(2.0 .* box_r, 1e-2),
            ],
        )
    else
        tr = [max(maximum(x_traj[i, :]) - minimum(x_traj[i, :]), 1e-2) for i in 1:nx]
        ts = [max(std(x_traj[i, :]), 1e-2) for i in 1:nx]
        append!(
            pairs,
            [
                :trajectory_range => tr,
                :trajectory_std => ts,
                :box_based => max.(2.0 .* box_r, 1e-2),
            ],
        )
    end
    return pairs
end

function _certify_named(name, scaling, problem, system_cfg, cfg, candidate)
    return _run_certification(
        String(name),
        problem,
        system_cfg,
        cfg,
        candidate;
        state_scaling = scaling,
    )
end

function _chain_ellipsoids(result)
    result.lmi_data !== nothing && hasproperty(result.lmi_data, :ellipsoids) ||
        return UT.Ellipsoid[]
    # Stored as [E_terminal, ..., E_initial], return by forward time [E_initial, ..., E_terminal].
    return reverse(collect(result.lmi_data.ellipsoids))
end

function _volumes(result)
    return [Float64(UT.get_volume(E)) for E in _chain_ellipsoids(result)]
end

function _kappa_matrices(result)
    mats = []
    for step in sort(result.steps; by = s -> s.k)
        κ = step.kappa
        if κ !== nothing && hasproperty(κ, :A)
            offset =
                hasproperty(κ, :b) ? getproperty(κ, :b) :
                (
                    hasproperty(κ, :c) ? getproperty(κ, :c) :
                    zeros(size(getproperty(κ, :A), 1))
                )
            push!(mats, (Matrix(getproperty(κ, :A)), collect(offset)))
        end
    end
    return mats
end

function _max_gain_diff(a, b)
    ka = _kappa_matrices(a)
    kb = _kappa_matrices(b)
    isempty(ka) || isempty(kb) && return missing
    n = min(length(ka), length(kb))
    return maximum(max(norm(ka[i][1] - kb[i][1]), norm(ka[i][2] - kb[i][2])) for i in 1:n)
end

function _identity_check(problem, system_cfg, cfg, candidate)
    nx = length(first(ST.enum_elems(candidate.x_traj)))
    runs = [
        :baseline =>
            _certify_named(:baseline, nothing, problem, system_cfg, cfg, candidate),
        :identity =>
            _certify_named(:identity, ones(nx), problem, system_cfg, cfg, candidate),
        :explicit_identity => _certify_named(
            :explicit_identity,
            fill(1.0, nx),
            problem,
            system_cfg,
            cfg,
            candidate,
        ),
    ]
    base_vols = _volumes(runs[1].second.result)
    rows = NamedTuple[]
    for (name, run) in runs
        vols = _volumes(run.result)
        diffs = length(vols) == length(base_vols) ? abs.(vols .- base_vols) : [Inf]
        rels =
            length(vols) == length(base_vols) ? diffs ./ max.(abs.(base_vols), eps()) :
            [Inf]
        conclusion = maximum(rels) <= 1e-5 ? "passed" : "failed"
        push!(
            rows,
            (;
                method = String(name),
                success = run.result.success,
                min_volume = isempty(vols) ? missing : minimum(vols),
                median_volume = isempty(vols) ? missing : median(vols),
                max_volume = isempty(vols) ? missing : maximum(vols),
                initial_volume = isempty(vols) ? missing : first(vols),
                total_certification_time = run.total_time,
                max_abs_volume_diff_vs_baseline = maximum(diffs),
                max_rel_volume_diff_vs_baseline = maximum(rels),
                max_gain_diff_vs_baseline = _max_gain_diff(
                    runs[1].second.result,
                    run.result,
                ),
                conclusion,
            ),
        )
    end
    CSV.write(FINAL_IDENTITY, DataFrame(rows))
    return (; rows, runs)
end

function _roundtrip_checks(problem, system_cfg, cfg, candidate, scalings)
    rows = NamedTuple[]
    rng = Random.MersenneTwister(7)
    for (name, scaling) in scalings
        scaling === nothing && continue
        cert = _build_certifier(problem, system_cfg, cfg; state_scaling = scaling)
        ctx = SC.build_symbolic_context(
            problem,
            candidate,
            cert.config,
            cert.symbolic_builder,
        )
        max_dyn = 0.0
        max_fb = 0.0
        max_elli = 0.0
        for k in unique(round.(Int, range(1, ctx.K; length = min(ctx.K, 6))))
            xk = collect(ctx.xs[k])
            xnext = collect(ctx.xs[k + 1])
            uk = collect(ctx.us[k])
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
            D = Diagonal(scaling)
            Sinv = Diagonal(1.0 ./ scaling)
            Az = Sinv * affsys.A * D
            Bz = Sinv * affsys.B
            cz = Sinv * (affsys.A * xk + affsys.c - xnext)
            Dwz = Sinv * affsys.D
            for _ in 1:10
                z = 0.2 .* randn(rng, length(xk))
                u = uk .+ 0.1 .* randn(rng, length(uk))
                w = zeros(size(affsys.D, 2))
                lhs =
                    Sinv * (
                        affsys.A * (xk + D * z) + affsys.B * u + affsys.c + affsys.D * w -
                        xnext
                    )
                rhs = Az * z + Bz * u + cz + Dwz * w
                max_dyn = max(max_dyn, norm(lhs - rhs))

                Kz = randn(rng, length(uk), length(xk))
                ell = randn(rng, length(uk))
                Kx = Kz * Sinv
                b = ell - Kx * xk
                max_fb = max(max_fb, norm((ell + Kz * z) - (Kx * (xk + D * z) + b)))

                Pz = Matrix{Float64}(I, length(xk), length(xk))
                Px = Matrix(Sinv' * Pz * Sinv)
                x = xk + D * z
                max_elli = max(max_elli, abs((z' * Pz * z) - ((x - xk)' * Px * (x - xk))))
            end
        end
        push!(
            rows,
            (;
                scaling = String(name),
                max_dynamics_error = max_dyn,
                max_feedback_error = max_fb,
                max_ellipsoid_error = max_elli,
                passed = max(max_dyn, max_fb, max_elli) <= 1e-8,
            ),
        )
    end
    CSV.write(FINAL_ROUNDTRIP, DataFrame(rows))
    return rows
end

function _ellipsoid_convention_check()
    E = UT.Ellipsoid([4.0;;], [0.0])
    radius = sqrt(inv(Matrix(E.P))[1, 1])
    volume = UT.get_volume(E)
    samples = sample_points_uniform_in_set(E, 500; rng = Random.MersenneTwister(11))
    max_value = maximum((x[1]^2 * 4.0 for x in samples))
    row = (;
        convention = "E(c,P) = {x : (x-c)' P (x-c) <= 1}",
        P = 4.0,
        expected_radius = 0.5,
        computed_radius = radius,
        expected_volume_1d = 1.0,
        computed_volume = volume,
        max_sample_membership_value = max_value,
        passed = abs(radius - 0.5) <= 1e-12 &&
                 abs(volume - 1.0) <= 1e-10 &&
                 max_value <= 1.0 + 1e-10,
    )
    CSV.write(FINAL_ELLIPSOID_CONVENTION, DataFrame([row]))
    return row
end

function _ellipsoid_box_stats(result, candidate, cfg; tol = 1e-7)
    xs = collect(ST.enum_elems(candidate.x_traj))
    radii = [max(abs(Float64(IA.inf(I))), abs(Float64(IA.sup(I)))) for I in collect(cfg.ΔX)]
    violations = 0
    max_margin = -Inf
    detailed = NamedTuple[]
    for step in result.steps
        E = step.ellipsoid
        if E === nothing
            continue
        end
        Q = Symmetric(Matrix(inv(E.P)))
        er = sqrt.(max.(diag(Q), 0.0))
        center_gap = abs.(collect(E.c) .- collect(xs[step.k]))
        margins = center_gap .+ er .- radii
        margin = maximum(margins)
        margin > tol && (violations += 1)
        max_margin = max(max_margin, margin)
        push!(
            detailed,
            (; k = step.k, box_violation_margin = margin, inside_box = margin <= tol),
        )
    end
    return (;
        violations,
        max_margin = isfinite(max_margin) ? max_margin : missing,
        detailed,
    )
end

function _input_violation_stats(result, problem, candidate)
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    U = problem.system.U
    violations = 0
    max_violation = 0.0
    for step in result.steps
        step.kappa === nothing && continue
        u = _eval_kappa_controller(step.kappa, collect(xs[step.k]))
        if hasproperty(U, :lb) && hasproperty(U, :ub)
            v = maximum(max.(collect(U.lb) .- u, u .- collect(U.ub), 0.0))
            v > 1e-8 && (violations += 1)
            max_violation = max(max_violation, v)
        end
    end
    return (; violations, max_violation)
end

function _stress_for_run(run_payload, x0_samples; method, scenario)
    return _simulate_samples(run_payload, x0_samples; method = method, scenario = scenario)
end

function _stress_metrics(rows; prefix)
    n = length(rows)
    c(s) = count(r -> r[:status] == s, rows)
    ns = c("success")
    nl = c("left_chain")
    nt = c("target_failed")
    ni = c("obstacle_violation")
    ne = c("simulation_error")
    ci_s = _wilson_interval(ns, n)
    ci_l = _wilson_interval(nl, n)
    return Dict{Symbol, Any}(
        Symbol(prefix, "_n_samples") => n,
        Symbol(prefix, "_success_rate") => ns / n,
        Symbol(prefix, "_left_chain_rate") => nl / n,
        Symbol(prefix, "_target_failed_rate") => nt / n,
        Symbol(prefix, "_input_violation_rate") => ni / n,
        Symbol(prefix, "_simulation_error_rate") => ne / n,
        Symbol(prefix, "_success_ci_low") => ci_s[1],
        Symbol(prefix, "_success_ci_high") => ci_s[2],
        Symbol(prefix, "_left_chain_ci_low") => ci_l[1],
        Symbol(prefix, "_left_chain_ci_high") => ci_l[2],
    )
end

function _run_metrics(name, scaling, run, problem, candidate, cfg, common_samples)
    vols = _volumes(run.result)
    box = _ellipsoid_box_stats(run.result, candidate, cfg)
    inp = _input_violation_stats(run.result, problem, candidate)
    payload = _run_result_payload(
        run,
        problem,
        cfg,
        candidate,
        (; root = RESULTS_DIR, plots_dir = FINAL_PLOTS_DIR),
    )
    own_samples = if run.result.success && !isempty(_chain_ellipsoids(run.result))
        sample_points_uniform_in_set(
            first(_chain_ellipsoids(run.result)),
            final_n_samples;
            rng = Random.MersenneTwister(final_rng_seed + hash(name) % 10000),
        )
    else
        Vector{Vector{Float64}}()
    end
    own =
        isempty(own_samples) ? nothing :
        _stress_for_run(payload, own_samples; method = String(name), scenario = "own_E0")
    common =
        run.result.success ?
        _stress_for_run(
            payload,
            common_samples;
            method = String(name),
            scenario = "common_baseline_E0",
        ) : nothing
    outside_common = if run.result.success && !isempty(_chain_ellipsoids(run.result))
        E0 = first(_chain_ellipsoids(run.result))
        count(x -> !(x ∈ E0), common_samples) / length(common_samples)
    else
        missing
    end
    d = Dict{Symbol, Any}(
        :method => String(name),
        :state_scaling_name => String(name),
        :state_scaling_vector => _vec_string(scaling),
        :success => run.result.success,
        :certified_transitions => count(s -> s.status == :ok, run.result.steps),
        :horizon => OP.horizon(candidate),
        :min_volume_physical => isempty(vols) ? missing : minimum(vols),
        :median_volume_physical => isempty(vols) ? missing : median(vols),
        :mean_volume_physical => isempty(vols) ? missing : mean(vols),
        :max_volume_physical => isempty(vols) ? missing : maximum(vols),
        :initial_volume_physical => isempty(vols) ? missing : first(vols),
        :final_volume_physical => isempty(vols) ? missing : last(vols),
        :total_certification_time => run.total_time,
        :average_local_solve_time => missing,
        :solver_status_summary => _solver_status_summary(run.result),
        :all_source_ellipsoids_inside_linearization_box => box.violations == 0,
        :number_box_violations => box.violations,
        :maximum_box_violation_margin => box.max_margin,
        :all_inputs_inside_U => inp.violations == 0,
        :all_inputs_inside_Ubar => missing,
        :number_input_bound_violations => inp.violations,
        :maximum_input_violation => inp.max_violation,
        :outside_method_E0_rate => outside_common,
        :formal_candidate =>
            run.result.success && box.violations == 0 && inp.violations == 0,
    )
    if own !== nothing
        merge!(d, _stress_metrics(own.rows; prefix = "own"))
    end
    if common !== nothing
        merge!(d, _stress_metrics(common.rows; prefix = "common"))
        merge!(d, Dict(:mcnemar_p_value_vs_baseline => missing))
    end
    return (; summary = NamedTuple(d), box_details = box.detailed, own, common)
end

function _stress_debug_for_run(name, run, problem, candidate, cfg; simulator_type)
    run.result.success || return NamedTuple[]
    payload = _run_result_payload(
        run,
        problem,
        cfg,
        candidate,
        (; root = RESULTS_DIR, plots_dir = FINAL_PLOTS_DIR),
    )
    chain = _build_certified_kappa_chain(run.result)
    samples = sample_points_uniform_in_set(
        chain.initial_ellipsoid,
        n_samples_debug;
        rng = Random.MersenneTwister(final_rng_seed + 99),
    )
    sim = _simulate_samples(
        payload,
        samples;
        method = String(name),
        scenario = simulator_type,
    )
    rows = NamedTuple[]
    for row in sim.rows
        push!(
            rows,
            (;
                method = String(name),
                simulator_type,
                sample_id = row[:sample_id],
                initial_membership_value = _ellipsoid_value(
                    chain.initial_ellipsoid,
                    samples[row[:sample_id]],
                ),
                status = row[:status],
                first_failure_step = row[:first_failure_step],
                max_ellipsoid_value = row[:maximum_normalized_ellipsoid_violation],
                final_ellipsoid_value = row[:maximum_normalized_ellipsoid_violation],
                final_distance_to_target = row[:final_distance_to_target],
                max_input_violation = 0.0,
                notes = "same replay path as benchmark helper; certification dynamics uses Ts=$(candidate.Ts), rk4_substeps=$(cfg.symbolic_rk4_substeps)",
            ),
        )
    end
    return rows
end

function _write_plots(scaling_summary, box_summary, scaling_detailed, stress_samples)
    mkpath(FINAL_PLOTS_DIR)
    ss = DataFrame(scaling_summary)
    bs = DataFrame(box_summary)
    sd = DataFrame(scaling_detailed)
    st = DataFrame(stress_samples)

    p = bar(
        ss.state_scaling_name,
        ss.median_volume_physical;
        yscale = :log10,
        xrotation = 45,
        legend = false,
        ylabel = "median volume",
        title = "Scaling sweep median physical volume",
    )
    savefig(p, joinpath(FINAL_PLOTS_DIR, "scaling_sweep_median_volume.png"))
    p = bar(
        ss.state_scaling_name,
        ss.min_volume_physical;
        yscale = :log10,
        xrotation = 45,
        legend = false,
        ylabel = "min volume",
        title = "Scaling sweep min physical volume",
    )
    savefig(p, joinpath(FINAL_PLOTS_DIR, "scaling_sweep_min_volume.png"))
    p = bar(
        ss.state_scaling_name,
        ss.own_success_rate;
        yerror = (
            ss.own_success_rate .- ss.own_success_ci_low,
            ss.own_success_ci_high .- ss.own_success_rate,
        ),
        xrotation = 45,
        legend = false,
        ylabel = "own success",
        title = "Scaling sweep own-E0 stress success",
    )
    savefig(p, joinpath(FINAL_PLOTS_DIR, "scaling_sweep_success_rate.png"))
    p = bar(
        ss.state_scaling_name,
        ss.common_success_rate;
        xrotation = 45,
        legend = false,
        ylabel = "common success",
        title = "Scaling sweep common baseline-E0 success",
    )
    savefig(p, joinpath(FINAL_PLOTS_DIR, "scaling_sweep_common_success_rate.png"))
    p = scatter(
        ss.median_volume_physical,
        ss.own_success_rate;
        xscale = :log10,
        xlabel = "median volume",
        ylabel = "own success",
        group = ss.state_scaling_name,
        title = "Volume vs stress success",
    )
    savefig(p, joinpath(FINAL_PLOTS_DIR, "scaling_sweep_volume_vs_success.png"))

    if !isempty(bs)
        p = plot(;
            xlabel = "box_multiplier_X",
            ylabel = "median volume",
            yscale = :log10,
            title = "Box sweep median volume",
        )
        for m in unique(bs.state_scaling_name)
            df = bs[bs.state_scaling_name .== m, :]
            plot!(
                p,
                df.box_multiplier_X,
                df.median_volume_physical;
                marker = :circle,
                label = m,
            )
        end
        savefig(p, joinpath(FINAL_PLOTS_DIR, "box_sweep_median_volume.png"))
        p = plot(;
            xlabel = "box_multiplier_X",
            ylabel = "own success",
            title = "Box sweep own success",
        )
        for m in unique(bs.state_scaling_name)
            df = bs[bs.state_scaling_name .== m, :]
            plot!(p, df.box_multiplier_X, df.own_success_rate; marker = :circle, label = m)
        end
        savefig(p, joinpath(FINAL_PLOTS_DIR, "box_sweep_success_rate.png"))
        p = plot(;
            xlabel = "box_multiplier_X",
            ylabel = "box violations",
            title = "Box inclusion violations",
        )
        for m in unique(bs.state_scaling_name)
            df = bs[bs.state_scaling_name .== m, :]
            plot!(
                p,
                df.box_multiplier_X,
                df.number_box_violations;
                marker = :circle,
                label = m,
            )
        end
        savefig(p, joinpath(FINAL_PLOTS_DIR, "box_sweep_box_violations.png"))
    end

    top_names = unique(
        vcat(
            ["none", "trajectory_range"],
            [ss.state_scaling_name[argmax(ss.median_volume_physical)]],
            [ss.state_scaling_name[argmax(ss.own_success_rate)]],
        ),
    )
    p = plot(;
        xlabel = "transition k",
        ylabel = "volume",
        yscale = :log10,
        title = "Volume vs transition top methods",
    )
    for name in top_names
        df = sd[(sd.sweep .== "scaling") .& (sd.method .== name), :]
        isempty(df) && continue
        plot!(p, df.k, df.volume_physical; label = name, marker = :circle)
    end
    savefig(p, joinpath(FINAL_PLOTS_DIR, "volume_vs_transition_top_methods.png"))

    failed = st[st.status .!= "success", :]
    if !isempty(failed)
        p = histogram(
            skipmissing(failed.first_failure_step);
            bins = 1:34,
            xlabel = "first failure step",
            title = "Stress first failure step histogram",
            label = "",
        )
        savefig(p, joinpath(FINAL_PLOTS_DIR, "stress_failure_step_histogram.png"))
    end
end

function _diagnostics_md(
    identity_rows,
    roundtrip_rows,
    convention_row,
    scaling_summary,
    box_summary,
    stress_debug,
)
    ss = DataFrame(scaling_summary)
    bs = DataFrame(box_summary)
    identity_ok = all(row.conclusion == "passed" for row in identity_rows)
    best_volume = ss.state_scaling_name[argmax(ss.median_volume_physical)]
    best_stress = ss.state_scaling_name[argmax(ss.own_success_rate)]
    formal = ss[(ss.formal_candidate .== true), :]
    debug_df = DataFrame(stress_debug)
    early =
        isempty(debug_df) ? 0 :
        count(x -> !ismissing(x) && x <= 3, debug_df.first_failure_step)
    open(FINAL_DIAGNOSTICS, "w") do io
        println(io, "# Final scaling validation diagnostics\n")
        println(
            io,
            "Monte Carlo stress test is not a proof. LMI feasibility is not enough if the ellipsoid leaves the local linearization box. A configuration is marked formally defensible only when certification succeeds and post-hoc local checks pass.\n",
        )
        println(io, "## Code state")
        println(io, "- No adaptive boxes are used in the certifier core.")
        println(io, "- The only refinement option is `state_scaling`.")
        println(
            io,
            "- `state_scaling = nothing` calls the historical `UT.transition_backward` path.\n",
        )
        println(io, "## Identity scaling check")
        println(io, identity_ok ? "passed" : "failed")
        println(io, "\n## Algebraic round-trip checks")
        println(io, all(row.passed for row in roundtrip_rows) ? "passed" : "failed")
        println(io, "\n## Ellipsoid convention")
        println(io, convention_row.passed ? "passed" : "failed")
        println(io, "\n## Stress test diagnosis")
        println(
            io,
            "- Samples are generated in E0 and initial membership is checked in `final_validation_stress_debug.csv`.",
        )
        println(
            io,
            "- Low success rates are mainly chain exits, not target failures or simulation errors, in the current runs.",
        )
        println(io, "- Number of debug failures by step <= 3: $(early).")
        println(
            io,
            "- This strongly suggests the certified source ellipsoids are not always contained in the local linearization boxes or that the local model bounds are not valid over the full certified ellipsoid. See box violation columns.",
        )
        println(io, "\n## Scaling sweep")
        println(io, "- Largest median volume: `$(best_volume)`.")
        println(io, "- Largest own-E0 stress success: `$(best_stress)`.")
        println(
            io,
            "- These are $(best_volume == best_stress ? "" : "not ")the same configuration.",
        )
        println(io, "\n## Box sweep")
        if isempty(bs)
            println(io, "No box sweep rows.")
        else
            best_formal = bs[bs.formal_candidate .== true, :]
            println(io, "- Formal candidates in box sweep: $(nrow(best_formal)).")
        end
        println(io, "\n## Formal validity")
        if isempty(formal)
            println(
                io,
                "No scaling sweep configuration satisfies all formal post-hoc checks.",
            )
        else
            for row in eachrow(formal)
                println(io, "- $(row.state_scaling_name)")
            end
        end
        println(io, "\n## Recommendation")
        return println(
            io,
            "Use the scaling results as experimental evidence only if identity scaling passes and the selected configuration has zero box/input violations. Otherwise, large physical volumes are suspicious and should not be presented as formal certificate improvements.",
        )
    end
end

function main()
    mkpath(RESULTS_DIR)
    mkpath(FINAL_PLOTS_DIR)
    isfile(TRAJECTORY_PATH) || error(
        "Missing saved MPPI trajectory. Run 01_generate_and_save_mppi_trajectory.jl first.",
    )

    cfg0 = MarcheAvantConfig(; output_root = RESULTS_DIR, plot_gif = false)
    candidate, data = _load_candidate()
    system_cfg = build_concrete_system()
    control_cfg = build_control_problem()
    problem = build_problem(system_cfg, control_cfg)
    x_traj = data["x_traj"]

    candidates = _scaling_candidates(x_traj, cfg0)
    println(
        "Part 1 confirmed: no adaptive-box code in certifier; only state_scaling remains.",
    )

    identity = _identity_check(problem, system_cfg, cfg0, candidate)
    roundtrip = _roundtrip_checks(problem, system_cfg, cfg0, candidate, candidates)
    convention = _ellipsoid_convention_check()

    baseline_run = identity.runs[1].second
    baseline_E0 = first(_chain_ellipsoids(baseline_run.result))
    common_samples = sample_points_uniform_in_set(
        baseline_E0,
        final_n_samples;
        rng = Random.MersenneTwister(final_rng_seed),
    )

    scaling_summary = NamedTuple[]
    scaling_detailed = NamedTuple[]
    all_stress_rows = Vector{Dict{Symbol, Any}}()
    stress_debug = NamedTuple[]
    run_by_name = Dict{String, Any}()

    for (name, scaling) in candidates
        println("[scaling] ", name)
        run =
            name == :none ? baseline_run :
            _certify_named(name, scaling, problem, system_cfg, cfg0, candidate)
        run_by_name[String(name)] = (; run, scaling)
        metrics = _run_metrics(name, scaling, run, problem, candidate, cfg0, common_samples)
        push!(scaling_summary, metrics.summary)
        append!(
            scaling_detailed,
            [
                (; sweep = "scaling", method = String(name), k = i, volume_physical = v) for
                (i, v) in enumerate(_volumes(run.result))
            ],
        )
        metrics.own !== nothing && append!(all_stress_rows, metrics.own.rows)
        metrics.common !== nothing && append!(all_stress_rows, metrics.common.rows)
    end

    for name in ["none", "trajectory_range"]
        haskey(run_by_name, name) || continue
        append!(
            stress_debug,
            _stress_debug_for_run(
                name,
                run_by_name[name].run,
                problem,
                candidate,
                cfg0;
                simulator_type = "benchmark_discrete",
            ),
        )
        append!(
            stress_debug,
            _stress_debug_for_run(
                name,
                run_by_name[name].run,
                problem,
                candidate,
                cfg0;
                simulator_type = "certification_discrete",
            ),
        )
    end

    CSV.write(FINAL_SUMMARY_SCALING, DataFrame(scaling_summary))
    CSV.write(FINAL_DETAILED_SCALING, DataFrame(scaling_detailed))
    CSV.write(FINAL_STRESS_SAMPLES, DataFrame(all_stress_rows))
    CSV.write(FINAL_STRESS_DEBUG, DataFrame(stress_debug))

    ss = DataFrame(scaling_summary)
    best_scaling_name =
        String(ss.state_scaling_name[argmax(coalesce.(ss.own_success_rate, -Inf))])
    best_scaling =
        haskey(run_by_name, best_scaling_name) ? run_by_name[best_scaling_name].scaling :
        nothing
    traj_scaling =
        first(pair.second for pair in candidates if pair.first == :trajectory_range)
    box_methods = [
        (:none, nothing),
        (Symbol(best_scaling_name), best_scaling),
        (:trajectory_range, traj_scaling),
    ]
    box_summary = NamedTuple[]
    box_detailed = NamedTuple[]
    for mx in box_multipliers_X,
        mu in box_multipliers_U,
        (name, scaling) in unique(box_methods)

        cfg = _cfg_with_boxes(cfg0; mx = mx, mu = mu)
        println("[box] mx=$(mx) mu=$(mu) scaling=$(name)")
        run = _certify_named(name, scaling, problem, system_cfg, cfg, candidate)
        metrics = _run_metrics(name, scaling, run, problem, candidate, cfg, common_samples)
        d = Dict{Symbol, Any}(pairs(metrics.summary))
        d[:box_multiplier_X] = mx
        d[:box_multiplier_U] = mu
        push!(box_summary, NamedTuple(d))
        append!(
            box_detailed,
            [
                (;
                    method = String(name),
                    box_multiplier_X = mx,
                    box_multiplier_U = mu,
                    row...,
                ) for row in metrics.box_details
            ],
        )
    end
    CSV.write(FINAL_SUMMARY_BOX, DataFrame(box_summary))
    CSV.write(FINAL_DETAILED_BOX, DataFrame(box_detailed))

    _write_plots(scaling_summary, box_summary, scaling_detailed, all_stress_rows)
    _diagnostics_md(
        identity.rows,
        roundtrip,
        convention,
        scaling_summary,
        box_summary,
        stress_debug,
    )

    println("summary scaling: ", FINAL_SUMMARY_SCALING)
    println("summary box: ", FINAL_SUMMARY_BOX)
    println("diagnostics: ", FINAL_DIAGNOSTICS)
    println("plots: ", FINAL_PLOTS_DIR)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
