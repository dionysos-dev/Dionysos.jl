using CSV
using DataFrames
using Plots
using Printf
using Statistics

ir_csv_transform(_col, val) = val === nothing ? missing : val
ir_csvwrite(path, table) = CSV.write(path, table; transform = ir_csv_transform)

function ir_cfg_tuple(cfg)
    names = fieldnames(typeof(cfg))
    return NamedTuple{names}(getfield(cfg, n) for n in names)
end

function ir_with_cfg(cfg; kwargs...)
    return typeof(cfg)(; merge(ir_cfg_tuple(cfg), NamedTuple(kwargs))...)
end

function ir_step_field(step, name::Symbol, default = missing)
    return hasproperty(step.summary, name) ? getproperty(step.summary, name) : default
end

function ir_step_value(step, name::Symbol, default = missing)
    value = ir_step_field(step, name, default)
    return value === nothing ? missing : value
end

function ir_principal_radii(E)
    vals = LA.eigvals(LA.Symmetric(Matrix{Float64}(E.P)))
    return sort(1.0 ./ sqrt.(max.(Float64.(vals), eps(Float64))); rev = true)
end

function ir_ellipsoid_chain(cert)
    rows = NamedTuple[]
    cert === nothing && return rows
    for step in cert.steps
        step.status == :ok || continue
        step.ellipsoid === nothing && continue
        push!(rows, (; k = step.k, E = step.ellipsoid))
    end
    return sort(rows; by = r -> r.k)
end

function ir_safe_ratio_margin(numer, denom)
    numer === nothing && return missing
    denom === nothing && return missing
    numer === missing && return missing
    denom === missing && return missing
    n = Float64.(collect(numer))
    d = Float64.(collect(denom))
    isempty(n) && return missing
    length(n) == length(d) || return missing
    return minimum([
        abs(b) <= eps(Float64) ? (a >= -eps(Float64) ? Inf : -Inf) : a / b for
        (a, b) in zip(n, d)
    ])
end

function ir_step_rows(label, alpha, cert)
    rows = NamedTuple[]
    cert === nothing && return rows
    for step in sort(cert.steps; by = s -> s.k)
        E = step.ellipsoid
        volume = E === nothing ? missing : Float64(UT.get_volume(E))
        radii = E === nothing ? Float64[] : ir_principal_radii(E)
        push!(
            rows,
            (;
                label,
                alpha,
                k = step.k,
                status = string(step.status),
                solver_status = string(ir_step_value(step, :solver_status, step.status)),
                volume,
                log_volume = volume === missing ? missing : log(volume),
                min_principal_radius = isempty(radii) ? missing : minimum(radii),
                min_X_containment_margin = ir_safe_ratio_margin(
                    ir_step_value(step, :Xbar_radius),
                    ir_step_value(step, :required_X_radius),
                ),
                min_U_containment_margin = ir_safe_ratio_margin(
                    ir_step_value(step, :Ubar_radius),
                    ir_step_value(step, :required_U_radius),
                ),
            ),
        )
    end
    return rows
end

function ir_control_vectors(candidate)
    candidate === nothing && return Vector{Vector{Float64}}()
    return [Float64.(collect(u)) for u in ST.enum_elems(candidate.u_traj)]
end

function ir_state_vectors(candidate)
    candidate === nothing && return Vector{Vector{Float64}}()
    return [Float64.(collect(x)) for x in ST.enum_elems(candidate.x_traj)]
end

function ir_nominal_reaches_target(problem, cfg, candidate; wrap_state)
    candidate === nothing && return false
    return any(x -> (wrap_state(x) ∈ problem.target_set), ST.enum_elems(candidate.x_traj))
end

function ir_input_reserve_metrics(
    problem,
    cfg,
    candidate;
    input_bounds_fn,
    planning_bounds_fn,
    tol::Float64 = 1.0e-8,
)
    us = ir_control_vectors(candidate)
    umin_cert, umax_cert = input_bounds_fn(problem, cfg)
    umin_plan, umax_plan = planning_bounds_fn(problem, cfg)
    cert_margins = [min(u[1] - umin_cert, umax_cert - u[1]) for u in us]
    plan_margins = [min(u[1] - umin_plan, umax_plan - u[1]) for u in us]
    return (;
        min_nominal_input_margin = isempty(cert_margins) ? missing : minimum(cert_margins),
        mean_nominal_input_margin = isempty(cert_margins) ? missing : mean(cert_margins),
        median_nominal_input_margin = isempty(cert_margins) ? missing :
                                      median(cert_margins),
        min_abs_distance_to_U_plan_boundary = isempty(plan_margins) ? missing :
                                              minimum(abs.(plan_margins)),
        number_of_steps_at_or_near_U_cert_boundary = count(m -> m <= tol, cert_margins),
        number_of_steps_at_or_near_U_plan_boundary = count(
            m -> abs(m) <= tol,
            plan_margins,
        ),
        umin_cert,
        umax_cert,
        umin_plan,
        umax_plan,
        input_boundary_tolerance = tol,
    )
end

function ir_certification_metrics(cert, candidate)
    chain = ir_ellipsoid_chain(cert)
    volumes = [Float64(UT.get_volume(r.E)) for r in chain]
    min_radii = [minimum(ir_principal_radii(r.E)) for r in chain]
    steps = ir_step_rows("", NaN, cert)
    xms = [
        r.min_X_containment_margin for
        r in steps if r.min_X_containment_margin !== missing &&
        isfinite(Float64(r.min_X_containment_margin))
    ]
    ums = [
        r.min_U_containment_margin for
        r in steps if r.min_U_containment_margin !== missing &&
        isfinite(Float64(r.min_U_containment_margin))
    ]
    return (;
        certification_complete = cert.success,
        certified_transitions = count(step -> step.status == :ok, cert.steps),
        total_transitions = candidate === nothing ? missing : OP.horizon(candidate),
        first_failed_backward_step = cert.failed_k === nothing ? missing : cert.failed_k,
        min_true_volume = isempty(volumes) ? missing : minimum(volumes),
        median_true_volume = isempty(volumes) ? missing : median(volumes),
        max_true_volume = isempty(volumes) ? missing : maximum(volumes),
        min_principal_radius = isempty(min_radii) ? missing : minimum(min_radii),
        median_min_principal_radius = isempty(min_radii) ? missing : median(min_radii),
        min_X_containment_margin = isempty(xms) ? missing : minimum(xms),
        min_U_containment_margin = isempty(ums) ? missing : minimum(ums),
        strict_box_containment = cert.success &&
                                 !isempty(xms) &&
                                 !isempty(ums) &&
                                 minimum(xms) >= 1.0 - 1.0e-8 &&
                                 minimum(ums) >= 1.0 - 1.0e-8,
    )
end

function ir_rollout_success_rate(stat_result)
    stat_result === nothing && return missing
    summary = stat_result.summary
    return hasproperty(summary, :final_target_rate) ? summary.final_target_rate : missing
end

function ir_failure_category(row)
    row.nominal_reaches_target == true || return "MPPI failed to reach target"
    row.certification_complete == true || return "certification incomplete"
    row.rollout_success_rate === missing && return "rollout validation not run"
    row.rollout_success_rate == 1.0 || return "statistical rollout validation failed"
    row.strict_box_containment == true || return "strict box containment violated"
    return "success"
end

function ir_save_step_csv(path, rows)
    mkpath(dirname(path))
    return ir_csvwrite(path, DataFrame(rows))
end

function ir_plot_campaign_diagnostics(results_dir, summary_df, all_step_rows; benchmark)
    plots_dir = joinpath(results_dir, "plots")
    mkpath(plots_dir)

    if nrow(summary_df) > 0
        fig_margin = plot(;
            xlabel = "alpha",
            ylabel = "margin to physical U_cert",
            title = "$(benchmark): nominal input reserve",
            legend = :topright,
            size = (720, 420),
        )
        plot!(
            fig_margin,
            summary_df.alpha,
            summary_df.min_nominal_input_margin;
            marker = :circle,
            label = "min",
        )
        plot!(
            fig_margin,
            summary_df.alpha,
            summary_df.mean_nominal_input_margin;
            marker = :diamond,
            label = "mean",
        )
        hline!(
            fig_margin,
            [0.0];
            color = :black,
            linestyle = :dash,
            label = "U_cert boundary",
        )
        savefig(fig_margin, joinpath(plots_dir, "$(benchmark)_input_margin_summary.pdf"))

        fig_trans = bar(
            string.(summary_df.alpha),
            summary_df.certified_transitions;
            xlabel = "alpha",
            ylabel = "certified transitions",
            title = "$(benchmark): certified transitions",
            legend = false,
            size = (720, 420),
        )
        savefig(
            fig_trans,
            joinpath(plots_dir, "$(benchmark)_certified_transitions_summary.pdf"),
        )
    end

    step_df = DataFrame(all_step_rows)
    if nrow(step_df) > 0
        fig_v = plot(;
            xlabel = "backward step k",
            ylabel = "true ellipsoid volume",
            yscale = :log10,
            title = "$(benchmark): true-scale volumes",
            legend = :topright,
            size = (720, 420),
        )
        fig_r = plot(;
            xlabel = "backward step k",
            ylabel = "minimum principal radius",
            yscale = :log10,
            title = "$(benchmark): true-scale radii",
            legend = :topright,
            size = (720, 420),
        )
        for alpha in unique(step_df.alpha)
            sub = step_df[step_df.alpha .== alpha, :]
            okv = sub[.!ismissing.(sub.volume), :]
            isempty(okv.k) ||
                plot!(fig_v, okv.k, okv.volume; marker = :circle, label = "alpha=$(alpha)")
            okr = sub[.!ismissing.(sub.min_principal_radius), :]
            isempty(okr.k) || plot!(
                fig_r,
                okr.k,
                okr.min_principal_radius;
                marker = :diamond,
                label = "alpha=$(alpha)",
            )
        end
        savefig(fig_v, joinpath(plots_dir, "$(benchmark)_true_volume_vs_step.pdf"))
        savefig(fig_r, joinpath(plots_dir, "$(benchmark)_true_min_radius_vs_step.pdf"))
    end
end

function ir_write_report(path, summary_df; benchmark)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# MPPI input reserve campaign: $(benchmark)\n")
        println(
            io,
            "This campaign tests whether small certified ellipsoids are primarily caused by nominal input saturation. The MPPI planner was constrained to use a strict subset of the physical input bounds, while the certifier kept the full admissible input set. Therefore, any improvement in certification can be interpreted as the effect of leaving input authority for the local feedback controller.\n",
        )
        if nrow(summary_df) == 0
            println(io, "No completed rows were produced.\n")
            return
        end
        baseline_rows = summary_df[summary_df.alpha .== 1.0, :]
        baseline = nrow(baseline_rows) == 0 ? nothing : baseline_rows[1, :]
        best_complete = summary_df[summary_df.certification_complete .== true, :]
        println(io, "## Summary")
        alpha_list = join(string.(summary_df.alpha), ", ")
        println(io, "- alpha values tested: $(alpha_list)")
        println(
            io,
            "- complete certified chains: $(sum(summary_df.certification_complete .== true)) / $(nrow(summary_df))",
        )
        println(
            io,
            "- nominal target successes: $(sum(summary_df.nominal_reaches_target .== true)) / $(nrow(summary_df))",
        )
        println(
            io,
            "- strict box-containment successes: $(sum(summary_df.strict_box_containment .== true)) / $(nrow(summary_df))\n",
        )

        if baseline !== nothing
            println(io, "## Baseline alpha = 1.0")
            println(
                io,
                "- nominal input margin min/mean/median: $(baseline.min_nominal_input_margin), $(baseline.mean_nominal_input_margin), $(baseline.median_nominal_input_margin)",
            )
            println(io, "- certification complete: $(baseline.certification_complete)")
            println(
                io,
                "- certified transitions: $(baseline.certified_transitions) / $(baseline.total_transitions)",
            )
            println(io, "- min true volume: $(baseline.min_true_volume)")
            println(io, "- min principal radius: $(baseline.min_principal_radius)")
            println(io, "- strict box containment: $(baseline.strict_box_containment)\n")
        end

        improved_margin =
            baseline === nothing ? false :
            any(
                row ->
                    row.alpha < 1.0 &&
                    row.min_nominal_input_margin !== missing &&
                    Float64(row.min_nominal_input_margin) >
                    Float64(baseline.min_nominal_input_margin),
                eachrow(summary_df),
            )
        complete_improved =
            baseline === nothing ? false :
            any(
                row ->
                    row.alpha < 1.0 &&
                    row.certification_complete == true &&
                    row.min_true_volume !== missing &&
                    baseline.min_true_volume !== missing &&
                    Float64(row.min_true_volume) > Float64(baseline.min_true_volume),
                eachrow(summary_df),
            )
        target_failures = sum(summary_df.nominal_reaches_target .== false)
        cert_failures = sum(summary_df.certification_complete .== false)

        println(io, "## Interpretation")
        println(
            io,
            "1. Did restricting MPPI input range increase nominal input margin? $(improved_margin).",
        )
        println(
            io,
            "2. Did restricted nominal trajectories still reach the target? $(target_failures == 0).",
        )
        println(
            io,
            "3. Did certification improve with a complete chain and larger true ellipsoids? $(complete_improved).",
        )
        println(
            io,
            "4. Did any chain become strictly sound? $(any(summary_df.strict_box_containment .== true)).",
        )
        println(
            io,
            "5. Incomplete certifications: $(cert_failures). Nominal target failures: $(target_failures).\n",
        )

        if complete_improved
            println(
                io,
                "Restricting the MPPI input range increased the nominal input margin and produced a complete certified chain with larger true-scale ellipsoids. This supports the hypothesis that planning with input reserve improves local certifiability.",
            )
        elseif target_failures > 0 && improved_margin
            println(
                io,
                "Restricting the MPPI input range increased the nominal input margin but prevented at least one nominal trajectory from reaching the target. This reveals a reachability/certifiability trade-off: the maneuver requires a large fraction of the available control authority, which makes robust local certification difficult.",
            )
        elseif cert_failures > 0
            println(
                io,
                "Restricting the MPPI input range did not produce a complete improved certificate in this bounded sweep. This suggests that input saturation is not the only bottleneck; nonlinear remainder bounds, box containment, anisotropy, or the local LMI formulation also limit the certified ellipsoid sizes.",
            )
        else
            println(
                io,
                "The input-reserve sweep did not show a clear improvement over alpha = 1.0 in the tested configurations.",
            )
        end
        return println(
            io,
            "\nAll ellipsoid volumes and radii reported here are true-scale diagnostics from the certifier output. No artificial plotting inflation is used.",
        )
    end
    return path
end
