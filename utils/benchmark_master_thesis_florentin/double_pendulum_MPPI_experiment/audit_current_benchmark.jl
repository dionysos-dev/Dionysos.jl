include(joinpath(@__DIR__, "audit_certification_boxes.jl"))

import Dates
import LinearAlgebra as LA
import Printf
import Statistics

const CURRENT_AUDIT_DIR = joinpath(@__DIR__, "outputs", "audit", "current_benchmark")

_diag_has(x, name::Symbol) = hasproperty(x, name)
_diag_get(x, name::Symbol, default = nothing) =
    _diag_has(x, name) ? getproperty(x, name) : default

function _diag_sigdigits(x; digits::Int = 5)
    x === nothing && return "nothing"
    x isa Real || return string(x)
    !isfinite(Float64(x)) && return string(x)
    return Printf.@sprintf("%.*g", digits, Float64(x))
end

function _diag_vec(v; digits::Int = 5)
    v === nothing && return "nothing"
    return "[" * join((_diag_sigdigits(x; digits) for x in collect(v)), ", ") * "]"
end

function _diag_csv_cell(x)
    s = x isa AbstractVector ? _diag_vec(x; digits = 8) : string(x)
    return occursin(',', s) || occursin('"', s) || occursin('\n', s) ?
           "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

function _diag_write_csv(path::AbstractString, rows)
    isempty(rows) && return path
    mkpath(dirname(path))
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for row in rows
            println(io, join((_diag_csv_cell(getproperty(row, n)) for n in names), ","))
        end
    end
    return path
end

function _diag_logvolume(E)
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

function _diag_axis_radii(E)
    E === nothing && return nothing
    try
        P = Matrix{Float64}(E.P)
        Q = inv(LA.Symmetric(P))
        return sqrt.(max.(0.0, LA.diag(Q)))
    catch
        return nothing
    end
end

function _diag_principal_radii(E)
    E === nothing && return nothing
    try
        vals = LA.eigvals(LA.Symmetric(Matrix{Float64}(E.P)))
        return sort(1.0 ./ sqrt.(max.(vals, eps(Float64))); rev = true)
    catch
        return nothing
    end
end

function _diag_max_ratio(required, box)
    required === nothing && return NaN
    box === nothing && return NaN
    r = Float64.(collect(required))
    b = Float64.(collect(box))
    isempty(r) && return NaN
    return maximum(r ./ max.(b, eps(Float64)))
end

function _diag_l_parts(L, nx::Int)
    L === nothing && return nothing, nothing
    Lv = Float64.(collect(L))
    return Lv[1:min(nx, length(Lv))], Lv[(min(nx, length(Lv)) + 1):end]
end

function _diag_interval_box_radii(box, n::Int)
    lows = Float64.(inf.(box))
    highs = Float64.(sup.(box))
    length(lows) == n || return fill(NaN, n)
    return 0.5 .* (highs .- lows)
end

function _diag_phase_label(depth, logvol, max_lx, max_ltail, max_rx_ratio, max_ru_ratio)
    if !isfinite(logvol)
        return "failed"
    elseif depth <= 3
        return "terminal-large"
    elseif max(max_rx_ratio, max_ru_ratio) > 1.0
        return "box-inconsistent"
    elseif max_lx > 1.0e-1 || max_ltail > 1.0e-1
        return "large-remainder"
    elseif logvol < -25.0
        return "tiny-ellipsoid"
    else
        return "regular"
    end
end

function _diag_rows(res, K::Int, nx::Int)
    rows = NamedTuple[]
    for (attempt, rec) in enumerate(res.steps)
        info = rec.summary
        E = rec.ellipsoid
        axis = _diag_axis_radii(E)
        principal = _diag_principal_radii(E)
        logvol = _diag_logvolume(E)
        Lx, Ltail = _diag_l_parts(_diag_get(info, :L), nx)
        Xbar = _diag_get(info, :Xbar_radius)
        Ubar = _diag_get(info, :Ubar_radius)
        reqX = _diag_get(info, :required_X_radius)
        reqU = _diag_get(info, :required_U_radius)
        max_rx_ratio = _diag_max_ratio(reqX, Xbar)
        max_ru_ratio = _diag_max_ratio(reqU, Ubar)
        max_lx = Lx === nothing || isempty(Lx) ? NaN : maximum(abs.(Lx))
        max_ltail = Ltail === nothing || isempty(Ltail) ? 0.0 : maximum(abs.(Ltail))
        depth = K - rec.k + 1
        push!(
            rows,
            (;
                attempt,
                k = rec.k,
                depth_from_terminal = depth,
                status = rec.status,
                adaptive_status = _diag_get(info, :adaptive_box_status),
                adaptive_iters = _diag_get(info, :adaptive_box_iters),
                selected_scale = _diag_get(info, :selected_scale),
                cost = rec.cost === nothing ? NaN : Float64(rec.cost),
                logvolume = logvol,
                volume = _audit_ellipsoid_volume(E),
                axis_radii = axis,
                min_axis_radius = axis === nothing ? NaN : minimum(axis),
                max_axis_radius = axis === nothing ? NaN : maximum(axis),
                principal_radii = principal,
                min_principal_radius = principal === nothing ? NaN : minimum(principal),
                max_principal_radius = principal === nothing ? NaN : maximum(principal),
                L_state = Lx,
                L_extra = Ltail,
                max_L_state = max_lx,
                max_L_extra = max_ltail,
                Xbar_radius = Xbar,
                Ubar_radius = Ubar,
                required_X_radius = reqX,
                required_U_radius = reqU,
                max_required_X_over_box = max_rx_ratio,
                max_required_U_over_box = max_ru_ratio,
                candidate_statuses = _diag_get(info, :candidate_statuses, Symbol[]),
                candidate_scales = _diag_get(info, :candidate_scales, Float64[]),
                phase = _diag_phase_label(
                    depth,
                    logvol,
                    max_lx,
                    max_ltail,
                    max_rx_ratio,
                    max_ru_ratio,
                ),
            ),
        )
    end
    return rows
end

function _diag_group_ranges(rows)
    isempty(rows) && return NamedTuple[]
    groups = NamedTuple[]
    start_i = 1
    current = rows[1].phase
    for i in 2:length(rows)
        if rows[i].phase != current
            chunk = rows[start_i:(i - 1)]
            push!(
                groups,
                (;
                    phase = current,
                    k_range = "$(first(chunk).k)-$(last(chunk).k)",
                    depth_range = "$(first(chunk).depth_from_terminal)-$(last(chunk).depth_from_terminal)",
                    count = length(chunk),
                    min_logvolume = minimum(r -> r.logvolume, chunk),
                    max_logvolume = maximum(r -> r.logvolume, chunk),
                    max_L_state = maximum(
                        r -> isfinite(r.max_L_state) ? r.max_L_state : -Inf,
                        chunk,
                    ),
                    max_required_X_over_box = maximum(
                        r ->
                            isfinite(r.max_required_X_over_box) ?
                            r.max_required_X_over_box : -Inf,
                        chunk,
                    ),
                    max_required_U_over_box = maximum(
                        r ->
                            isfinite(r.max_required_U_over_box) ?
                            r.max_required_U_over_box : -Inf,
                        chunk,
                    ),
                ),
            )
            start_i = i
            current = rows[i].phase
        end
    end
    chunk = rows[start_i:end]
    push!(
        groups,
        (;
            phase = current,
            k_range = "$(first(chunk).k)-$(last(chunk).k)",
            depth_range = "$(first(chunk).depth_from_terminal)-$(last(chunk).depth_from_terminal)",
            count = length(chunk),
            min_logvolume = minimum(r -> r.logvolume, chunk),
            max_logvolume = maximum(r -> r.logvolume, chunk),
            max_L_state = maximum(
                r -> isfinite(r.max_L_state) ? r.max_L_state : -Inf,
                chunk,
            ),
            max_required_X_over_box = maximum(
                r -> isfinite(r.max_required_X_over_box) ? r.max_required_X_over_box : -Inf,
                chunk,
            ),
            max_required_U_over_box = maximum(
                r -> isfinite(r.max_required_U_over_box) ? r.max_required_U_over_box : -Inf,
                chunk,
            ),
        ),
    )
    return groups
end

function _diag_top(rows, field::Symbol; n::Int = 5, finite_only::Bool = true)
    vals = [(i, getproperty(rows[i], field)) for i in eachindex(rows)]
    if finite_only
        vals = filter(x -> x[2] isa Real && isfinite(Float64(x[2])), vals)
    end
    sort!(vals; by = x -> Float64(x[2]), rev = true)
    return rows[[i for (i, _) in first(vals, min(n, length(vals)))]]
end

function _diag_write_report(path, cfg, candidate, res, rows, groups, state_scaling)
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    ok_rows = filter(r -> r.status == :ok, rows)
    failed_rows = filter(r -> r.status != :ok, rows)
    min_row = isempty(ok_rows) ? nothing : ok_rows[argmin([r.logvolume for r in ok_rows])]
    max_req_x =
        isempty(rows) ? NaN :
        maximum(
            r -> isfinite(r.max_required_X_over_box) ? r.max_required_X_over_box : -Inf,
            rows,
        )
    max_req_u =
        isempty(rows) ? NaN :
        maximum(
            r -> isfinite(r.max_required_U_over_box) ? r.max_required_U_over_box : -Inf,
            rows,
        )
    max_l_state =
        isempty(rows) ? NaN :
        maximum(r -> isfinite(r.max_L_state) ? r.max_L_state : -Inf, rows)
    max_l_extra =
        isempty(rows) ? NaN :
        maximum(r -> isfinite(r.max_L_extra) ? r.max_L_extra : -Inf, rows)

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Double Pendulum Current Benchmark Diagnostic")
        println(io)
        println(io, "- generated_at: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "- success: ", res.success)
        println(io, "- failed_k: ", res.failed_k)
        println(io, "- attempted_steps: ", length(rows))
        println(io, "- ok_steps: ", length(ok_rows))
        println(io, "- candidate_horizon: ", length(us))
        println(io, "- candidate_terminal: ", _diag_vec(xs[end]))
        println(io, "- candidate_last_u: ", _diag_vec(us[end]))
        println(io, "- state_scaling_mode: ", cfg.state_scaling_mode)
        println(
            io,
            "- state_scaling_diag: ",
            state_scaling === nothing ? "nothing" : _diag_vec(LA.diag(state_scaling)),
        )
        println(io, "- adaptive_linearization_boxes: ", cfg.adaptive_linearization_boxes)
        println(
            io,
            "- fixed_DeltaX_radius: ",
            _diag_vec(_diag_interval_box_radii(cfg.ΔX, 4)),
        )
        println(
            io,
            "- fixed_DeltaU_radius: ",
            _diag_vec(_diag_interval_box_radii(cfg.ΔU, 1)),
        )
        println(io, "- adaptive_DeltaX_initial: ", _diag_vec(cfg.ΔX_initial))
        println(io, "- adaptive_DeltaU_initial: ", _diag_vec(cfg.ΔU_initial))
        println(io, "- max_required_X_over_box: ", _diag_sigdigits(max_req_x))
        println(io, "- max_required_U_over_box: ", _diag_sigdigits(max_req_u))
        println(io, "- max_L_state: ", _diag_sigdigits(max_l_state))
        println(io, "- max_L_extra: ", _diag_sigdigits(max_l_extra))
        if min_row !== nothing
            println(
                io,
                "- smallest_ellipsoid_step: k=",
                min_row.k,
                " depth=",
                min_row.depth_from_terminal,
            )
            println(
                io,
                "- smallest_ellipsoid_logvolume: ",
                _diag_sigdigits(min_row.logvolume),
            )
            println(io, "- smallest_ellipsoid_axis_radii: ", _diag_vec(min_row.axis_radii))
            println(
                io,
                "- smallest_ellipsoid_principal_radii: ",
                _diag_vec(min_row.principal_radii),
            )
        end
        if !isempty(failed_rows)
            r = first(failed_rows)
            println(io, "- failure_status: ", r.status)
            println(io, "- failure_adaptive_status: ", r.adaptive_status)
            println(io, "- failure_Xbar_radius: ", _diag_vec(r.Xbar_radius))
            println(io, "- failure_Ubar_radius: ", _diag_vec(r.Ubar_radius))
            println(io, "- failure_L_state: ", _diag_vec(r.L_state))
            println(io, "- failure_L_extra: ", _diag_vec(r.L_extra))
            println(io, "- failure_required_X_radius: ", _diag_vec(r.required_X_radius))
            println(io, "- failure_required_U_radius: ", _diag_vec(r.required_U_radius))
        end

        println(io)
        println(io, "## Phase Summary")
        println(io)
        println(
            io,
            "| phase | k range | depth range | count | logvol min | logvol max | max Lx | max reqX/box | max reqU/box |",
        )
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for g in groups
            println(
                io,
                "| ",
                g.phase,
                " | ",
                g.k_range,
                " | ",
                g.depth_range,
                " | ",
                g.count,
                " | ",
                _diag_sigdigits(g.min_logvolume),
                " | ",
                _diag_sigdigits(g.max_logvolume),
                " | ",
                _diag_sigdigits(g.max_L_state),
                " | ",
                _diag_sigdigits(g.max_required_X_over_box),
                " | ",
                _diag_sigdigits(g.max_required_U_over_box),
                " |",
            )
        end

        println(io)
        println(io, "## Largest Remainder Bounds")
        println(io)
        println(io, "| k | depth | phase | max Lx | L state | L extra | logvol |")
        println(io, "|---:|---:|---|---:|---|---|---:|")
        for r in _diag_top(rows, :max_L_state; n = 8)
            println(
                io,
                "| ",
                r.k,
                " | ",
                r.depth_from_terminal,
                " | ",
                r.phase,
                " | ",
                _diag_sigdigits(r.max_L_state),
                " | ",
                _diag_vec(r.L_state),
                " | ",
                _diag_vec(r.L_extra),
                " | ",
                _diag_sigdigits(r.logvolume),
                " |",
            )
        end

        println(io)
        println(io, "## Largest Box Consistency Ratios")
        println(io)
        println(
            io,
            "| k | depth | phase | reqX/box | reqU/box | reqX | Xbar | reqU | Ubar |",
        )
        println(io, "|---:|---:|---|---:|---:|---|---|---|---|")
        for r in _diag_top(rows, :max_required_U_over_box; n = 8)
            println(
                io,
                "| ",
                r.k,
                " | ",
                r.depth_from_terminal,
                " | ",
                r.phase,
                " | ",
                _diag_sigdigits(r.max_required_X_over_box),
                " | ",
                _diag_sigdigits(r.max_required_U_over_box),
                " | ",
                _diag_vec(r.required_X_radius),
                " | ",
                _diag_vec(r.Xbar_radius),
                " | ",
                _diag_vec(r.required_U_radius),
                " | ",
                _diag_vec(r.Ubar_radius),
                " |",
            )
        end
    end
    return path
end

function run_current_double_pendulum_diagnostic(;
    cfg::DoublePendulumMPPIConfig = DoublePendulumMPPIConfig(;
        plot_gif = false,
        plot_mp4 = false,
        verbose = false,
        load_saved_trajectory = true,
    ),
    output_dir::AbstractString = CURRENT_AUDIT_DIR,
)
    mkpath(output_dir)
    raw_candidate = load_candidate_trajectory(cfg.saved_trajectory_path)
    candidate = default_prepare_for_certification(cfg)(raw_candidate)
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))

    problem = build_double_pendulum_problem(cfg; pendulum_module = DP)
    system_cfg = build_double_pendulum_system_cfg(cfg; pendulum_module = DP)
    control_cfg = build_double_pendulum_control_cfg(cfg; pendulum_module = DP)
    cert = build_double_pendulum_certifier(problem, system_cfg, control_cfg, cfg)
    SC.set_problem!(cert, problem)
    SC.set_trajectory!(cert, candidate)

    println("candidate_horizon=", length(us))
    println("candidate_terminal=", collect(Float64, xs[end]))
    println("candidate_last_u=", collect(Float64, us[end]))
    println(
        "state_scaling=",
        cfg.state_scaling_mode,
        " ",
        LA.diag(resolve_state_scaling(cfg, candidate)),
    )
    println("adaptive_linearization_boxes=", cfg.adaptive_linearization_boxes)
    println("running certification diagnostic...")
    t = @elapsed SC.certify!(cert)
    res = SC.get_result(cert)
    println(
        "done success=",
        res.success,
        " failed_k=",
        res.failed_k,
        " attempted_steps=",
        length(res.steps),
        " solve_time=",
        round(t; digits = 2),
    )

    rows = _diag_rows(res, length(us), length(xs[1]))
    groups = _diag_group_ranges(rows)
    csv_path = joinpath(output_dir, "current_config_step_diagnostics.csv")
    phase_path = joinpath(output_dir, "current_config_phase_summary.csv")
    report_path = joinpath(output_dir, "current_config_diagnostic.md")
    _diag_write_csv(csv_path, rows)
    _diag_write_csv(phase_path, groups)
    _diag_write_report(
        report_path,
        cfg,
        candidate,
        res,
        rows,
        groups,
        resolve_state_scaling(cfg, candidate),
    )

    println("step_csv=", csv_path)
    println("phase_csv=", phase_path)
    println("report=", report_path)
    return (;
        cfg,
        candidate,
        cert,
        result = res,
        rows,
        groups,
        csv_path,
        phase_path,
        report_path,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_current_double_pendulum_diagnostic()
end
