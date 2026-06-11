include(joinpath(@__DIR__, "run_simple_pendulum_mppi.jl"))

const AUDIT_ROOT = joinpath(@__DIR__, "outputs", "audit")

function _audit_vec_string(x)
    x === nothing && return ""
    return replace(string(collect(x)), "," => ";")
end

function _audit_csv_cell(x)
    s = string(x)
    return occursin(',', s) || occursin('"', s) || occursin('\n', s) ?
           "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

function _audit_step_metric(cert, name::Symbol; reducer = identity, default = NaN)
    vals = Float64[]
    for step in cert.steps
        step.status == :ok || continue
        val = _step_field(step, name, nothing)
        val === nothing && continue
        try
            y = reducer(val)
            y === nothing && continue
            push!(vals, Float64(y))
        catch
        end
    end
    isempty(vals) && return default
    return sum(vals) / length(vals)
end

function _audit_ellipsoid_volume(E)
    E === nothing && return NaN
    try
        return Float64(UT.get_volume(E))
    catch
        return NaN
    end
end

function _audit_initial_ellipsoid(cert)
    cert === nothing && return nothing
    cert.lmi_data === nothing && return nothing
    hasproperty(cert.lmi_data, :ellipsoids) || return nothing
    isempty(cert.lmi_data.ellipsoids) && return nothing
    return last(cert.lmi_data.ellipsoids)
end

function _audit_terminal_dists(run_result, cfg)
    terminal_data = terminal_inner_ellipsoid_data(run_result.problem, cfg)
    xs = collect(ST.enum_elems(run_result.nominal_candidate.x_traj))
    terminal_data === nothing && return (;
        final_dT2 = NaN,
        min_dT2 = NaN,
        hit_index = 0,
        endpoint = _audit_vec_string(xs[end]),
    )

    dists = [
        terminal_ellipsoidal_distance2(
            build_pendulum_wrap_state(cfg)(x),
            terminal_data.terminal_center,
            terminal_data.terminal_shape,
        ) for x in xs
    ]
    hit = findfirst(d -> d <= cfg.terminal_success_distance2 + 1.0e-8, dists)
    return (;
        final_dT2 = Float64(last(dists)),
        min_dT2 = Float64(minimum(dists)),
        hit_index = hit === nothing ? 0 : hit,
        endpoint = _audit_vec_string(xs[end]),
    )
end

function _audit_input_metrics(run_result)
    us = collect(ST.enum_elems(run_result.nominal_candidate.u_traj))
    isempty(us) && return (; max_abs_u = NaN, rms_u = NaN, total_variation_u = NaN)
    uvals = [Float64(u[1]) for u in us]
    max_abs_u = maximum(abs, uvals)
    rms_u = sqrt(sum(abs2, uvals) / length(uvals))
    total_variation_u =
        length(uvals) <= 1 ? 0.0 :
        sum(abs(uvals[i] - uvals[i - 1]) for i in 2:length(uvals))
    return (; max_abs_u, rms_u, total_variation_u)
end

function _audit_row(name::AbstractString, cfg, run_result)
    cert = run_result.result.certification
    terminal = _audit_terminal_dists(run_result, cfg)
    inputs = _audit_input_metrics(run_result)
    E0 = _audit_initial_ellipsoid(cert)
    stat = run_result.statistical
    summary = stat === nothing ? nothing : stat.summary
    active_scaling =
        run_result.solver.certifier isa SimplePendulumScalingCertifier ?
        active_state_scaling(run_result.solver.certifier) : certifier_state_scaling(cfg)

    return (;
        name,
        success = cert.success,
        failed_k = cert.failed_k === nothing ? "" : cert.failed_k,
        cert_steps = length(cert.steps),
        candidate_horizon = OP.horizon(run_result.nominal_candidate),
        terminal_hit_index = terminal.hit_index,
        final_dT2 = terminal.final_dT2,
        min_dT2 = terminal.min_dT2,
        endpoint = terminal.endpoint,
        initial_volume = _audit_ellipsoid_volume(E0),
        lmi_sec = run_result.timings.lmi_sec,
        mppi_sec = run_result.timings.mppi_sec,
        total_sec = run_result.timings.total_sec,
        max_abs_u = inputs.max_abs_u,
        rms_u = inputs.rms_u,
        total_variation_u = inputs.total_variation_u,
        stat_certified_success_rate = summary === nothing ? NaN :
                                      summary.certified_success_rate,
        stat_closed_loop_success_rate = summary === nothing ? NaN :
                                        summary.closed_loop_success_rate,
        stat_ellipsoid_exit_rate = summary === nothing ? NaN : summary.ellipsoid_exit_rate,
        avg_selected_logvolume = _audit_step_metric(cert, :selected_logvolume),
        avg_L = _audit_step_metric(cert, :L),
        avg_box_theta = _audit_step_metric(
            cert,
            :Xbar_radius;
            reducer = x -> length(x) >= 1 ? x[1] : nothing,
        ),
        avg_box_omega = _audit_step_metric(
            cert,
            :Xbar_radius;
            reducer = x -> length(x) >= 2 ? x[2] : nothing,
        ),
        avg_box_u = _audit_step_metric(
            cert,
            :Ubar_radius;
            reducer = x -> length(x) >= 1 ? x[1] : nothing,
        ),
        planning_input_scale = cfg.planning_input_scale,
        mppi_noise_u = cfg.mppi_noise_u,
        terminal_shrink = cfg.terminal_shrink,
        terminal_success_distance2 = cfg.terminal_success_distance2,
        state_scaling_mode = cfg.state_scaling_mode,
        state_scaling = _audit_vec_string(active_scaling),
        adaptive_linearization_boxes = cfg.adaptive_linearization_boxes,
    )
end

function _write_audit_csv(path::AbstractString, rows)
    isempty(rows) && return path
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for row in rows
            println(io, join((_audit_csv_cell(getproperty(row, n)) for n in names), ","))
        end
    end
    return path
end

function _best_successful(rows)
    ok = [r for r in rows if r.success]
    isempty(ok) && return nothing
    sort!(
        ok;
        by = r -> (
            r.initial_volume,
            r.stat_certified_success_rate,
            -r.total_variation_u,
            -r.total_sec,
        ),
        rev = true,
    )
    return first(ok)
end

function _write_audit_report(path::AbstractString, rows)
    best = _best_successful(rows)
    open(path, "w") do io
        println(io, "# Simple pendulum MPPI option audit")
        println(io)
        println(
            io,
            "Rows are generated with plots/animations disabled and a reduced statistical sample.",
        )
        println(io)
        if best !== nothing
            println(
                io,
                "Best successful row by initial ellipsoid volume: `",
                best.name,
                "`.",
            )
            println(io)
        end
        println(
            io,
            "| name | success | E0 volume | cert stat | min dT2 | hit k | mppi s | lmi s | scaling | input scale |",
        )
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---|---:|")
        for r in rows
            println(
                io,
                "| ",
                r.name,
                " | ",
                r.success,
                " | ",
                round(r.initial_volume; sigdigits = 5),
                " | ",
                round(r.stat_certified_success_rate; sigdigits = 4),
                " | ",
                round(r.min_dT2; sigdigits = 5),
                " | ",
                r.terminal_hit_index,
                " | ",
                round(r.mppi_sec; sigdigits = 4),
                " | ",
                round(r.lmi_sec; sigdigits = 4),
                " | ",
                r.state_scaling,
                " | ",
                r.planning_input_scale,
                " |",
            )
        end
    end
    return path
end

function run_simple_pendulum_mppi_option_audit()
    mkpath(AUDIT_ROOT)

    cases = [
        ("baseline_current", (;)),
        ("terminal_threshold_1", (; terminal_success_distance2 = 1.0)),
        (
            "terminal_full_ellipsoid",
            (; terminal_shrink = 1.0, terminal_success_distance2 = 1.0),
        ),
        ("input_scale_0_8", (; planning_input_scale = 0.8, mppi_noise_u = 0.6)),
        ("input_scale_0_6", (; planning_input_scale = 0.6, mppi_noise_u = 0.5)),
        ("input_scale_0_7", (; planning_input_scale = 0.7, mppi_noise_u = 0.55)),
        ("input_scale_0_65", (; planning_input_scale = 0.65, mppi_noise_u = 0.5)),
        ("input_scale_0_55", (; planning_input_scale = 0.55, mppi_noise_u = 0.45)),
        ("input_scale_0_5", (; planning_input_scale = 0.5, mppi_noise_u = 0.4)),
        (
            "scaling_targetish",
            (;
                state_scaling_matrix = Matrix{Float64}(LA.Diagonal([0.30, 1.0])),
                state_scaling = [0.30, 1.0],
            ),
        ),
        (
            "scaling_domain",
            (;
                state_scaling_matrix = Matrix{Float64}(LA.Diagonal([pi, 7.0])),
                state_scaling = [pi, 7.0],
            ),
        ),
        ("fixed_boxes", (; adaptive_linearization_boxes = false)),
        (
            "combined_terminal_input_scale",
            (;
                terminal_shrink = 1.0,
                terminal_success_distance2 = 1.0,
                planning_input_scale = 0.8,
                mppi_noise_u = 0.6,
                state_scaling_matrix = Matrix{Float64}(LA.Diagonal([pi, 7.0])),
                state_scaling = [pi, 7.0],
            ),
        ),
        (
            "input_0_6_terminal_threshold_1",
            (;
                planning_input_scale = 0.6,
                mppi_noise_u = 0.5,
                terminal_success_distance2 = 1.0,
            ),
        ),
        (
            "input_0_6_terminal_full",
            (;
                planning_input_scale = 0.6,
                mppi_noise_u = 0.5,
                terminal_shrink = 1.0,
                terminal_success_distance2 = 1.0,
            ),
        ),
        (
            "input_0_6_terminal_full_center_0",
            (;
                planning_input_scale = 0.6,
                mppi_noise_u = 0.5,
                terminal_shrink = 1.0,
                terminal_success_distance2 = 1.0,
                terminal_center_weight = 0.0,
            ),
        ),
        (
            "input_0_6_terminal_full_center_1e3",
            (;
                planning_input_scale = 0.6,
                mppi_noise_u = 0.5,
                terminal_shrink = 1.0,
                terminal_success_distance2 = 1.0,
                terminal_center_weight = 1.0e3,
            ),
        ),
        (
            "input_0_6_terminal_full_center_1e4",
            (;
                planning_input_scale = 0.6,
                mppi_noise_u = 0.5,
                terminal_shrink = 1.0,
                terminal_success_distance2 = 1.0,
                terminal_center_weight = 1.0e4,
            ),
        ),
        (
            "input_0_6_scaling_targetish",
            (;
                planning_input_scale = 0.6,
                mppi_noise_u = 0.5,
                state_scaling_matrix = Matrix{Float64}(LA.Diagonal([0.30, 1.0])),
                state_scaling = [0.30, 1.0],
            ),
        ),
        (
            "input_0_6_scaling_domain",
            (;
                planning_input_scale = 0.6,
                mppi_noise_u = 0.5,
                state_scaling_matrix = Matrix{Float64}(LA.Diagonal([pi, 7.0])),
                state_scaling = [pi, 7.0],
            ),
        ),
        (
            "input_0_6_lambda_1e_4",
            (; planning_input_scale = 0.6, mppi_noise_u = 0.5, λ = 1.0e-4),
        ),
        (
            "input_0_6_lambda_1e_3",
            (; planning_input_scale = 0.6, mppi_noise_u = 0.5, λ = 1.0e-3),
        ),
    ]

    rows = NamedTuple[]
    for (name, params) in cases
        println("running audit case: ", name)
        cfg_params = merge(
            (;
                output_root = joinpath(AUDIT_ROOT, name),
                plot_gif = false,
                plot_mp4 = false,
                adaptive_box_verbose = false,
                kappa_statistical_samples = 80,
                kappa_statistical_seed = 11,
                rng_seed = 11,
            ),
            params,
        )
        cfg = SimplePendulumMPPIConfig(; cfg_params...)
        try
            run_result = redirect_stdout(devnull) do
                return main(
                    cfg;
                    save_outputs = false,
                    run_statistical = true,
                    artifact_prefix = name,
                )
            end
            push!(rows, _audit_row(name, cfg, run_result))
        catch err
            @warn "audit case failed" name exception = (err, catch_backtrace())
            push!(
                rows,
                (;
                    name,
                    success = false,
                    failed_k = "exception",
                    cert_steps = 0,
                    candidate_horizon = 0,
                    terminal_hit_index = 0,
                    final_dT2 = NaN,
                    min_dT2 = NaN,
                    endpoint = "",
                    initial_volume = NaN,
                    lmi_sec = NaN,
                    mppi_sec = NaN,
                    total_sec = NaN,
                    max_abs_u = NaN,
                    rms_u = NaN,
                    total_variation_u = NaN,
                    stat_certified_success_rate = NaN,
                    stat_closed_loop_success_rate = NaN,
                    stat_ellipsoid_exit_rate = NaN,
                    avg_selected_logvolume = NaN,
                    avg_L = NaN,
                    avg_box_theta = NaN,
                    avg_box_omega = NaN,
                    avg_box_u = NaN,
                    planning_input_scale = NaN,
                    mppi_noise_u = NaN,
                    terminal_shrink = NaN,
                    terminal_success_distance2 = NaN,
                    state_scaling_mode = :failed,
                    state_scaling = "",
                    adaptive_linearization_boxes = false,
                ),
            )
        end
    end

    csv_path = _write_audit_csv(
        joinpath(AUDIT_ROOT, "simple_pendulum_mppi_option_audit.csv"),
        rows,
    )
    report_path = _write_audit_report(
        joinpath(AUDIT_ROOT, "simple_pendulum_mppi_option_audit.md"),
        rows,
    )
    println("audit_csv = ", csv_path)
    println("audit_report = ", report_path)
    return (; rows, csv_path, report_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_simple_pendulum_mppi_option_audit()
end
