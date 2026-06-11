include(joinpath(@__DIR__, "run_double_pendulum_mppi.jl"))

import LinearAlgebra as LA

function _audit_suffix_candidate(cand, n::Integer)
    xs = collect(ST.enum_elems(cand.x_traj))
    us = collect(ST.enum_elems(cand.u_traj))
    K = length(us)
    n = min(Int(n), K)
    return OP.CandidateTrajectory(
        ST.Trajectory(xs[(K - n + 1):(K + 1)]),
        ST.Trajectory(us[(K - n + 1):K]);
        Ts = cand.Ts,
        source = :certification_audit_suffix,
        metadata = (; source = cand.source, suffix_steps = n),
    )
end

function _audit_interval_box(radii::AbstractVector{<:Real})
    return IA.IntervalBox((IA.interval(-Float64(r), Float64(r)) for r in radii)...)
end

function _audit_vec_string(v)
    v === nothing && return "nothing"
    return "[" * join(round.(Float64.(collect(v)); sigdigits = 4), ", ") * "]"
end

function _audit_csv_cell(x)
    s = string(x)
    return occursin(',', s) || occursin('"', s) || occursin('\n', s) ?
           "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

function _audit_append_csv_row!(path::AbstractString, row; header::Bool = false)
    mkpath(dirname(path))
    names = propertynames(row)
    open(path, header ? "w" : "a") do io
        if header
            println(io, join(string.(names), ","))
        end
        return println(io, join((_audit_csv_cell(getproperty(row, n)) for n in names), ","))
    end
    return path
end

function _audit_last_step_summary(res)
    isempty(res.steps) && return (;)
    rec = last(res.steps)
    info = rec.summary
    return (;
        k = rec.k,
        status = rec.status,
        adaptive_status = hasproperty(info, :adaptive_box_status) ?
                          info.adaptive_box_status : nothing,
        iters = hasproperty(info, :adaptive_box_iters) ? info.adaptive_box_iters : nothing,
        selected_scale = hasproperty(info, :selected_scale) ? info.selected_scale : nothing,
        Xbar_radius = hasproperty(info, :Xbar_radius) ? info.Xbar_radius : nothing,
        Ubar_radius = hasproperty(info, :Ubar_radius) ? info.Ubar_radius : nothing,
        required_X_radius = hasproperty(info, :required_X_radius) ? info.required_X_radius :
                            nothing,
        required_U_radius = hasproperty(info, :required_U_radius) ? info.required_U_radius :
                            nothing,
        candidate_statuses = hasproperty(info, :candidate_statuses) ?
                             info.candidate_statuses : Symbol[],
        candidate_scales = hasproperty(info, :candidate_scales) ? info.candidate_scales :
                           Float64[],
    )
end

function _audit_ellipsoid_volume(E)
    E === nothing && return NaN
    try
        return Float64(UT.get_volume(E))
    catch
        return NaN
    end
end

function _audit_ellipsoid_logvolume(E)
    E === nothing && return NaN
    try
        vol = _audit_ellipsoid_volume(E)
        isfinite(vol) && vol > 0.0 && return log(vol)
    catch
    end

    try
        P = Matrix{Float64}(E.P)
        F = LA.cholesky(LA.Symmetric(P))
        n = size(P, 1)
        unit_ball_logvol = (n / 2) * log(pi) - loggamma(n / 2 + 1)
        return Float64(unit_ball_logvol - sum(log, LA.diag(F.U)))
    catch
        return NaN
    end
end

function _audit_certified_ellipsoids(res)
    res === nothing && return Any[]
    res.lmi_data === nothing && return Any[]
    hasproperty(res.lmi_data, :ellipsoids) || return Any[]
    return collect(res.lmi_data.ellipsoids)
end

function _audit_volume_summary(res)
    ellipsoids = _audit_certified_ellipsoids(res)
    logs = [_audit_ellipsoid_logvolume(E) for E in ellipsoids]
    logs = [x for x in logs if isfinite(x)]
    ok_steps = count(step -> step.status == :ok, res.steps)

    if isempty(logs)
        return (;
            certified_steps = ok_steps,
            ellipsoid_count = 0,
            terminal_logvolume = NaN,
            deepest_logvolume = NaN,
            min_logvolume = NaN,
            mean_logvolume = NaN,
            max_logvolume = NaN,
        )
    end

    return (;
        certified_steps = ok_steps,
        ellipsoid_count = length(ellipsoids),
        terminal_logvolume = first(logs),
        deepest_logvolume = last(logs),
        min_logvolume = minimum(logs),
        mean_logvolume = sum(logs) / length(logs),
        max_logvolume = maximum(logs),
    )
end

function _audit_scaling_matrix(alpha::Real)
    a = Float64(alpha)
    return Matrix{Float64}(LA.Diagonal([pi^a, pi^a, 6.0^a, 6.0^a]))
end

function _audit_state_scaling_matrix(angle_scale::Real, velocity_scale::Real)
    return Matrix{Float64}(
        LA.Diagonal([
            Float64(angle_scale),
            Float64(angle_scale),
            Float64(velocity_scale),
            Float64(velocity_scale),
        ]),
    )
end

function _audit_reach_volume_score(row; reach_weight::Float64 = 100.0)
    logvol = isfinite(row.deepest_logvolume) ? row.deepest_logvolume : -Inf
    return reach_weight * Float64(row.certified_steps) + logvol
end

function _audit_variant_config(name::AbstractString)
    common = (;
        plot_gif = false,
        plot_mp4 = false,
        verbose = false,
        load_saved_trajectory = true,
        adaptive_box_verbose = false,
    )

    name == "default_config" && return DoublePendulumMPPIConfig(; common...)

    name == "default_more_iters" &&
        return DoublePendulumMPPIConfig(; common..., adaptive_box_max_iters = 12)

    name == "default_safety_103" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_box_safety = 1.03,
        adaptive_box_max_iters = 12,
    )

    name == "default_safety_105" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_box_safety = 1.05,
        adaptive_box_max_iters = 12,
    )

    name == "adaptive_current_fast" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_linearization_boxes = true,
        adaptive_box_max_iters = 1,
        adaptive_box_search_scales = [1.0],
        adaptive_box_keep_first_consistent = true,
    )

    name == "adaptive_small_fast" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_linearization_boxes = true,
        ΔX_initial = [0.02, 0.02, 0.02, 0.02],
        ΔU_initial = [0.02],
        adaptive_box_max_iters = 2,
        adaptive_box_search_scales = [1.0],
        adaptive_box_keep_first_consistent = true,
    )

    name == "adaptive_simple_style" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_linearization_boxes = true,
        ΔX_initial = [0.02, 0.02, 0.02, 0.02],
        ΔU_initial = [0.02],
        adaptive_box_safety = 1.01,
        adaptive_box_max_iters = 4,
        adaptive_box_search_scales = [0.7, 1.0, 1.5],
        adaptive_box_keep_first_consistent = false,
    )

    name == "adaptive_domain_scaling" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_linearization_boxes = true,
        state_scaling_mode = :matrix,
        state_scaling_matrix = Matrix{Float64}(LA.Diagonal([pi, pi, 6.0, 6.0])),
        ΔX_initial = [0.02, 0.02, 0.05, 0.05],
        ΔU_initial = [0.02],
        adaptive_box_safety = 1.01,
        adaptive_box_max_iters = 3,
        adaptive_box_search_scales = [1.0],
        adaptive_box_keep_first_consistent = true,
    )

    name == "fixed_005" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_linearization_boxes = false,
        ΔX = _audit_interval_box([0.05, 0.05, 0.05, 0.05]),
        ΔU = _audit_interval_box([0.05]),
    )

    name == "fixed_020" && return DoublePendulumMPPIConfig(;
        common...,
        adaptive_linearization_boxes = false,
        ΔX = _audit_interval_box([0.20, 0.20, 0.20, 0.20]),
        ΔU = _audit_interval_box([0.20]),
    )

    return error("Unknown audit variant $(name)")
end

function run_double_pendulum_scaling_volume_audit(;
    suffix_steps = [5],
    alphas = [-1.0, -0.5, 0.0, 0.5, 1.0],
    fixed_box_radius_x = [0.02, 0.02, 0.02, 0.02],
    fixed_box_radius_u = [0.02],
    use_fixed_boxes::Bool = true,
    adaptive_box_max_iters::Int = 12,
    adaptive_box_safety::Float64 = 1.03,
)
    base_cfg = DoublePendulumMPPIConfig(; plot_gif = false, plot_mp4 = false)
    candidate = default_prepare_for_certification(base_cfg)(
        load_candidate_trajectory(base_cfg.saved_trajectory_path),
    )
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    println("candidate_horizon=", length(us))
    println("candidate_terminal=", collect(Float64, xs[end]))
    println("candidate_last_u=", collect(Float64, us[end]))
    println("scaling_family=S(alpha)=diag([pi^alpha, pi^alpha, 6^alpha, 6^alpha])")
    println(
        "boxes=",
        use_fixed_boxes ? "fixed_unsound" : "adaptive_sound",
        " fixed_x=",
        fixed_box_radius_x,
        " fixed_u=",
        fixed_box_radius_u,
    )

    rows = NamedTuple[]
    for suffix in suffix_steps
        for alpha in alphas
            S = _audit_scaling_matrix(alpha)
            cfg_kwargs = (;
                plot_gif = false,
                plot_mp4 = false,
                verbose = false,
                load_saved_trajectory = true,
                state_scaling_mode = :matrix,
                state_scaling_matrix = S,
                adaptive_box_max_iters,
                adaptive_box_safety,
            )
            cfg = if use_fixed_boxes
                DoublePendulumMPPIConfig(;
                    cfg_kwargs...,
                    adaptive_linearization_boxes = false,
                    ΔX = _audit_interval_box(fixed_box_radius_x),
                    ΔU = _audit_interval_box(fixed_box_radius_u),
                )
            else
                DoublePendulumMPPIConfig(; cfg_kwargs...)
            end

            problem = build_double_pendulum_problem(cfg; pendulum_module = DP)
            system_cfg = build_double_pendulum_system_cfg(cfg; pendulum_module = DP)
            control_cfg = build_double_pendulum_control_cfg(cfg; pendulum_module = DP)
            sub = _audit_suffix_candidate(candidate, suffix)
            cert = build_double_pendulum_certifier(problem, system_cfg, control_cfg, cfg)
            SC.set_problem!(cert, problem)
            SC.set_trajectory!(cert, sub)

            print("alpha=", alpha, " scaling=", LA.diag(S), " suffix=", suffix, " ... ")
            flush(stdout)
            t = @elapsed SC.certify!(cert)
            res = SC.get_result(cert)
            vol = _audit_volume_summary(res)
            summary = _audit_last_step_summary(res)
            row = (;
                alpha = Float64(alpha),
                scaling = collect(Float64, LA.diag(S)),
                suffix = Int(suffix),
                success = res.success,
                failed_k = res.failed_k,
                solve_time = t,
                vol...,
                summary...,
            )
            push!(rows, row)
            println(
                "success=",
                row.success,
                " certified_steps=",
                row.certified_steps,
                " failed_k=",
                row.failed_k,
                " time=",
                round(row.solve_time; digits = 2),
                " deepest_logvol=",
                round(row.deepest_logvolume; digits = 4),
                " min_logvol=",
                round(row.min_logvolume; digits = 4),
                " status=",
                get(row, :status, nothing),
                " adaptive=",
                get(row, :adaptive_status, nothing),
            )
        end
    end

    comparable = filter(r -> isfinite(r.deepest_logvolume), rows)
    if !isempty(comparable)
        best_reach = maximum(r -> r.certified_steps, comparable)
        best_reach_rows = filter(r -> r.certified_steps == best_reach, comparable)
        best_volume =
            best_reach_rows[argmax([r.deepest_logvolume for r in best_reach_rows])]
        println(
            "best_among_max_reach alpha=",
            best_volume.alpha,
            " certified_steps=",
            best_volume.certified_steps,
            " success=",
            best_volume.success,
            " deepest_logvol=",
            round(best_volume.deepest_logvolume; digits = 4),
            " scaling=",
            best_volume.scaling,
        )
    end
    return rows
end

function run_double_pendulum_joint_scaling_box_audit(;
    suffix_steps = [10],
    state_angle_scales = [1 / pi, 1.0, sqrt(pi), pi],
    state_velocity_scales = [1 / 6.0, 1.0, sqrt(6.0), 6.0],
    input_scales = [0.5, 1.0, 2.0],
    state_box_scales = [0.5, 1.0, 2.0],
    input_box_scales = [0.5, 1.0, 2.0],
    base_box_radius_x = [0.02, 0.02, 0.02, 0.02],
    base_box_radius_u = [0.02],
    use_fixed_boxes::Bool = true,
    adaptive_box_max_iters::Int = 12,
    adaptive_box_safety::Float64 = 1.03,
    max_cases::Union{Nothing, Int} = nothing,
    results_path::String = joinpath(
        @__DIR__,
        "outputs",
        "audit",
        "double_pendulum_joint_scaling_box_audit.csv",
    ),
)
    base_cfg = DoublePendulumMPPIConfig(; plot_gif = false, plot_mp4 = false)
    candidate = default_prepare_for_certification(base_cfg)(
        load_candidate_trajectory(base_cfg.saved_trajectory_path),
    )
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    println("candidate_horizon=", length(us))
    println("candidate_terminal=", collect(Float64, xs[end]))
    println("candidate_last_u=", collect(Float64, us[end]))
    println("results_path=", results_path)
    println("boxes=", use_fixed_boxes ? "fixed_unsound" : "adaptive_sound")

    rows = NamedTuple[]
    wrote_header = false
    case_idx = 0
    for suffix in suffix_steps
        for angle_scale in state_angle_scales
            for velocity_scale in state_velocity_scales
                Sx = _audit_state_scaling_matrix(angle_scale, velocity_scale)
                for input_scale in input_scales
                    for xbox_scale in state_box_scales
                        for ubox_scale in input_box_scales
                            case_idx += 1
                            if max_cases !== nothing && case_idx > max_cases
                                return rows
                            end

                            x_box = Float64(xbox_scale) .* Float64.(base_box_radius_x)
                            u_box =
                                Float64(ubox_scale * input_scale) .*
                                Float64.(base_box_radius_u)
                            cfg_kwargs = (;
                                plot_gif = false,
                                plot_mp4 = false,
                                verbose = false,
                                load_saved_trajectory = true,
                                state_scaling_mode = :matrix,
                                state_scaling_matrix = Sx,
                                input_scaling_mode = :matrix,
                                input_scaling_matrix = Matrix{Float64}(
                                    LA.Diagonal([Float64(input_scale)]),
                                ),
                                maxδu = base_cfg.maxδu * Float64(input_scale)^2,
                                adaptive_box_max_iters,
                                adaptive_box_safety,
                            )
                            cfg = if use_fixed_boxes
                                DoublePendulumMPPIConfig(;
                                    cfg_kwargs...,
                                    adaptive_linearization_boxes = false,
                                    ΔX = _audit_interval_box(x_box),
                                    ΔU = _audit_interval_box(u_box),
                                )
                            else
                                DoublePendulumMPPIConfig(;
                                    cfg_kwargs...,
                                    adaptive_linearization_boxes = true,
                                    ΔX_initial = x_box,
                                    ΔU_initial = u_box,
                                )
                            end

                            problem =
                                build_double_pendulum_problem(cfg; pendulum_module = DP)
                            system_cfg =
                                build_double_pendulum_system_cfg(cfg; pendulum_module = DP)
                            control_cfg =
                                build_double_pendulum_control_cfg(cfg; pendulum_module = DP)
                            sub = _audit_suffix_candidate(candidate, suffix)
                            cert = build_double_pendulum_certifier(
                                problem,
                                system_cfg,
                                control_cfg,
                                cfg,
                            )
                            SC.set_problem!(cert, problem)
                            SC.set_trajectory!(cert, sub)

                            print(
                                "case=",
                                case_idx,
                                " suffix=",
                                suffix,
                                " sx=(",
                                angle_scale,
                                ",",
                                velocity_scale,
                                ") input_scale=",
                                input_scale,
                                " box=(",
                                xbox_scale,
                                ",",
                                ubox_scale,
                                ") ... ",
                            )
                            flush(stdout)

                            t = @elapsed SC.certify!(cert)
                            res = SC.get_result(cert)
                            vol = _audit_volume_summary(res)
                            summary = _audit_last_step_summary(res)
                            row = (;
                                case_idx,
                                suffix = Int(suffix),
                                angle_scale = Float64(angle_scale),
                                velocity_scale = Float64(velocity_scale),
                                input_scale = Float64(input_scale),
                                xbox_scale = Float64(xbox_scale),
                                ubox_scale = Float64(ubox_scale),
                                x_box = x_box,
                                u_box = u_box,
                                maxδu = cfg.maxδu,
                                success = res.success,
                                failed_k = res.failed_k,
                                solve_time = t,
                                vol...,
                                score = NaN,
                                summary...,
                            )
                            row = merge(row, (; score = _audit_reach_volume_score(row),))
                            push!(rows, row)
                            _audit_append_csv_row!(
                                results_path,
                                row;
                                header = !wrote_header,
                            )
                            wrote_header = true

                            println(
                                "success=",
                                row.success,
                                " certified_steps=",
                                row.certified_steps,
                                " failed_k=",
                                row.failed_k,
                                " deepest_logvol=",
                                round(row.deepest_logvolume; digits = 4),
                                " score=",
                                round(row.score; digits = 4),
                                " time=",
                                round(row.solve_time; digits = 2),
                                " status=",
                                get(row, :status, nothing),
                            )
                        end
                    end
                end
            end
        end
    end

    comparable = filter(r -> isfinite(r.deepest_logvolume), rows)
    if !isempty(comparable)
        best_reach = maximum(r -> r.certified_steps, comparable)
        best_reach_rows = filter(r -> r.certified_steps == best_reach, comparable)
        best_volume =
            best_reach_rows[argmax([r.deepest_logvolume for r in best_reach_rows])]
        best_score = comparable[argmax([r.score for r in comparable])]
        println(
            "best_volume_among_max_reach case=",
            best_volume.case_idx,
            " reach=",
            best_volume.certified_steps,
            " logvol=",
            round(best_volume.deepest_logvolume; digits = 4),
            " sx=(",
            best_volume.angle_scale,
            ",",
            best_volume.velocity_scale,
            ") input_scale=",
            best_volume.input_scale,
            " box=(",
            best_volume.xbox_scale,
            ",",
            best_volume.ubox_scale,
            ")",
        )
        println(
            "best_score case=",
            best_score.case_idx,
            " reach=",
            best_score.certified_steps,
            " logvol=",
            round(best_score.deepest_logvolume; digits = 4),
            " score=",
            round(best_score.score; digits = 4),
            " sx=(",
            best_score.angle_scale,
            ",",
            best_score.velocity_scale,
            ") input_scale=",
            best_score.input_scale,
            " box=(",
            best_score.xbox_scale,
            ",",
            best_score.ubox_scale,
            ")",
        )
    end
    return rows
end

function run_double_pendulum_certification_box_audit(;
    variants = [
        "adaptive_current_fast",
        "adaptive_small_fast",
        "adaptive_simple_style",
        "adaptive_domain_scaling",
        "fixed_005",
        "fixed_020",
    ],
    suffix_steps = [1],
)
    base_cfg = DoublePendulumMPPIConfig(; plot_gif = false, plot_mp4 = false)
    candidate = default_prepare_for_certification(base_cfg)(
        load_candidate_trajectory(base_cfg.saved_trajectory_path),
    )
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    println("candidate_horizon=", length(us))
    println("candidate_terminal=", collect(Float64, xs[end]))
    println("candidate_last_u=", collect(Float64, us[end]))

    rows = NamedTuple[]
    for variant in variants
        for suffix in suffix_steps
            cfg = _audit_variant_config(variant)
            problem = build_double_pendulum_problem(cfg; pendulum_module = DP)
            system_cfg = build_double_pendulum_system_cfg(cfg; pendulum_module = DP)
            control_cfg = build_double_pendulum_control_cfg(cfg; pendulum_module = DP)
            sub = _audit_suffix_candidate(candidate, suffix)
            cert = build_double_pendulum_certifier(problem, system_cfg, control_cfg, cfg)
            SC.set_problem!(cert, problem)
            SC.set_trajectory!(cert, sub)

            print("variant=", variant, " suffix=", suffix, " ... ")
            flush(stdout)
            t = @elapsed SC.certify!(cert)
            res = SC.get_result(cert)
            summary = _audit_last_step_summary(res)
            row = (;
                variant,
                suffix = Int(suffix),
                success = res.success,
                failed_k = res.failed_k,
                solve_time = t,
                summary...,
            )
            push!(rows, row)
            println(
                "success=",
                row.success,
                " failed_k=",
                row.failed_k,
                " time=",
                round(row.solve_time; digits = 2),
                " status=",
                get(row, :status, nothing),
                " adaptive=",
                get(row, :adaptive_status, nothing),
                " Xbar=",
                _audit_vec_string(get(row, :Xbar_radius, nothing)),
                " Ubar=",
                _audit_vec_string(get(row, :Ubar_radius, nothing)),
                " reqX=",
                _audit_vec_string(get(row, :required_X_radius, nothing)),
                " reqU=",
                _audit_vec_string(get(row, :required_U_radius, nothing)),
            )
        end
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    variants =
        haskey(ENV, "DP_CERT_AUDIT_VARIANTS") ? split(ENV["DP_CERT_AUDIT_VARIANTS"], ",") :
        ["adaptive_current_fast", "adaptive_small_fast"]
    suffix_steps =
        haskey(ENV, "DP_CERT_AUDIT_SUFFIXES") ?
        parse.(Int, split(ENV["DP_CERT_AUDIT_SUFFIXES"], ",")) : [1]
    run_double_pendulum_certification_box_audit(;
        variants = String.(variants),
        suffix_steps = suffix_steps,
    )
end
