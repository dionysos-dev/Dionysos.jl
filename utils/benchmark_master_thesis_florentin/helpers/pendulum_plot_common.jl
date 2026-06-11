function _pendulum_artifact_filename(
    basename::AbstractString,
    stem::AbstractString,
    ext::AbstractString = "pdf",
)
    name = isempty(basename) ? stem : "$(basename)_$(stem)"
    return "$(name).$(ext)"
end

function _prepare_periodic_state_plot_data(
    problem,
    candidate;
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims,
    periodic_periods,
    periodic_start,
    periodic_ellipsoid_copy_mode::Symbol = :all,
)
    xs = collect(ST.enum_elems(candidate.x_traj))
    if unwrap_angles
        xs = unwrap_periodic_state_list(xs, periodic_dims, periodic_periods)
    end
    if wrap_angles
        xs = wrap_periodic_state_list(xs, periodic_dims, periodic_periods, periodic_start)
    end

    raw_domain =
        (hasproperty(problem, :system) && hasproperty(problem.system, :X)) ?
        problem.system.X : nothing
    raw_initial_set = hasproperty(problem, :initial_set) ? problem.initial_set : nothing
    raw_target_set = hasproperty(problem, :target_set) ? problem.target_set : nothing

    domain = maybe_wrap_set(
        raw_domain,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )
    initial_set = maybe_wrap_set(
        raw_initial_set,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )
    target_set = maybe_wrap_set(
        raw_target_set,
        periodic_dims,
        periodic_periods,
        periodic_start;
        enabled = wrap_angles,
    )

    ellipsoids =
        (show_ellipsoids && cert_result !== nothing) ?
        extract_ellipsoids(cert_result; max_keep = nothing) : UT.Ellipsoid[]
    if !isempty(ellipsoids) && unwrap_angles
        ellipsoids = unwrap_ellipsoid_centers(ellipsoids, periodic_dims, periodic_periods)
    end
    if !isempty(ellipsoids) && wrap_angles
        ellipsoids = wrap_ellipsoid_centers(
            ellipsoids,
            periodic_dims,
            periodic_periods,
            periodic_start,
        )
        ellipsoids = add_periodic_ellipsoid_plot_copies(
            ellipsoids,
            periodic_dims,
            periodic_periods,
            periodic_start;
            enabled = true,
            mode = periodic_ellipsoid_copy_mode,
        )
    end

    us = collect(ST.enum_elems(candidate.u_traj))
    return (;
        xs,
        us,
        x_traj = ST.Trajectory(xs),
        raw_domain,
        raw_initial_set,
        raw_target_set,
        domain,
        initial_set,
        target_set,
        ellipsoids,
        wrap_angles,
        periodic_dims,
        periodic_periods,
        periodic_start,
        ts_x = candidate.Ts .* collect(0:(length(xs) - 1)),
        ts_u = candidate.Ts .* collect(0:(length(us) - 1)),
    )
end

function _plot_data_sets(data; use_raw_sets::Bool = false)
    if use_raw_sets
        return (;
            domain = data.raw_domain,
            initial_set = data.raw_initial_set,
            target_set = data.raw_target_set,
        )
    end
    return (;
        domain = data.domain,
        initial_set = data.initial_set,
        target_set = data.target_set,
    )
end

function _add_projected_ellipsoids!(fig, ellipsoids, dims; opacity::Float64 = 0.84)
    label_used = false
    for E in ellipsoids
        plot!(
            fig,
            project_ellipsoid_2d(E; dims = dims);
            color = THESIS_ELLIPSOID_COLOR,
            linecolor = PLOT_COLORS[:ellipsoid_edge],
            fillcolor = PLOT_COLORS[:ellipsoid],
            opacity = opacity,
            lw = 1.0,
            label = label_used ? "" : L"\mathrm{certified\ ellipsoids}",
        )
        label_used = true
    end
    return fig
end

function _projection_limits_from_data(
    data,
    dims;
    use_raw_sets::Bool = false,
    margin_ratio::Float64 = 0.08,
)
    sets = _plot_data_sets(data; use_raw_sets = use_raw_sets)
    bounds = Any[
        _trajectory_bounds_2d(data.xs, dims),
        _plot_bounds_2d(sets.initial_set, dims),
        _plot_bounds_2d(sets.target_set, dims),
    ]
    append!(bounds, [_projected_ellipsoid_bounds(E, dims) for E in data.ellipsoids])
    lims = _expand_bounds_2d(_merge_bounds_2d(bounds); margin_ratio = margin_ratio)

    if hasproperty(data, :wrap_angles) && data.wrap_angles
        xlims, ylims = lims.xlims, lims.ylims
        for (axis, dim) in enumerate(dims)
            lookup = _periodic_dim_lookup(dim, data.periodic_dims, data.periodic_periods)
            lookup === nothing && continue
            idx, period = lookup
            lo = Float64(data.periodic_start[idx])
            hi = lo + Float64(period)
            axis == 1 ? (xlims = (lo, hi)) : (ylims = (lo, hi))
        end
        return (; xlims, ylims)
    end

    return lims
end

function _build_state_space_projection(
    data;
    dims = (1, 2),
    title::AbstractString = "",
    use_raw_sets::Bool = false,
    ellipsoid_opacity::Float64 = 0.84,
    xlabel::AbstractString = "",
    ylabel::AbstractString = "",
    trajectory_label::AbstractString = "nominal trajectory",
    legend = :topright,
    size = (640, 460),
)
    sets = _plot_data_sets(data; use_raw_sets = use_raw_sets)
    lims = _projection_limits_from_data(data, dims; use_raw_sets = use_raw_sets)
    fig = plot(;
        thesis_plot_kwargs(; legend = legend, size = size)...,
        xlims = lims.xlims,
        ylims = lims.ylims,
        xlabel = xlabel,
        ylabel = ylabel,
    )
    plot_state_space_basic!(
        fig,
        sets.domain,
        sets.initial_set,
        sets.target_set,
        data.x_traj;
        dims = collect(dims),
        trajectory_label = trajectory_label,
        show_trajectory = false,
    )
    _plot_nominal_trajectory!(
        fig,
        data.x_traj;
        dims = collect(dims),
        label = trajectory_label,
    )
    _add_projected_ellipsoids!(fig, data.ellipsoids, dims; opacity = ellipsoid_opacity)
    return fig
end

function _build_phase_prefix_plot(
    data,
    k::Int;
    dims = (1, 2),
    title::AbstractString = "",
    use_raw_sets::Bool = false,
    xlabel::AbstractString = "",
    ylabel::AbstractString = "",
)
    fig = _build_state_space_projection(
        data;
        dims = dims,
        title = title,
        use_raw_sets = use_raw_sets,
        ellipsoid_opacity = 0.84,
        xlabel = xlabel,
        ylabel = ylabel,
        legend = false,
    )
    d1, d2 = dims
    prefix_1 = [x[d1] for x in data.xs[1:k]]
    prefix_2 = [x[d2] for x in data.xs[1:k]]
    plot!(
        fig,
        prefix_1,
        prefix_2;
        color = THESIS_TRAJ_COLOR,
        lw = 1.8,
        markershape = :circle,
        ms = 1.8,
        markerstrokewidth = 0,
        alpha = 0.95,
        label = "",
    )
    scatter!(
        fig,
        [data.xs[k][d1]],
        [data.xs[k][d2]];
        color = PLOT_COLORS[:target_edge],
        ms = 2.6,
        label = "",
    )
    return fig
end

function _add_input_bounds!(fig, problem)
    if hasproperty(problem.system, :U) && problem.system.U isa UT.HyperRectangle
        hline!(
            fig,
            [problem.system.U.lb[1], problem.system.U.ub[1]];
            color = PLOT_COLORS[:constraint],
            ls = :dash,
            label = "",
        )
    end
    return fig
end

function _build_control_time_plot(data, problem; title::AbstractString, k = nothing)
    fig = plot(
        data.ts_u,
        [u[1] for u in data.us];
        color = PLOT_COLORS[:nominal],
        lw = 1.4,
        markershape = :circle,
        ms = 1.8,
        markerstrokewidth = 0,
        label = L"u",
        xlabel = L"t\,[\mathrm{s}]",
        ylabel = L"u",
        thesis_plot_kwargs(; legend = :topright, size = (640, 360))...,
    )
    if k !== nothing && !isempty(data.us)
        uk_index = clamp(k, 1, length(data.us))
        scatter!(
            fig,
            [data.ts_u[uk_index]],
            [data.us[uk_index][1]];
            color = BENCH_ELLIPSOID_COLOR,
            ms = 5,
            label = "",
        )
    end
    return _add_input_bounds!(fig, problem)
end

function _save_state_time_plot!(
    output_path::AbstractString,
    data,
    problem,
    specs;
    title::AbstractString,
    layout,
    size,
)
    fig = plot(; layout = layout, thesis_plot_kwargs(; legend = :topright, size = size)...)
    target_lb = problem.target_set.lb
    target_ub = problem.target_set.ub

    for (i, spec) in enumerate(specs)
        values = [x[spec.index] for x in data.xs]
        plot!(
            fig[i],
            data.ts_x,
            values;
            color = spec.color,
            lw = 0.9,
            markershape = :circle,
            ms = 1.8,
            markerstrokewidth = 0,
            label = spec.label,
            xlabel = L"t\,[\mathrm{s}]",
            ylabel = spec.ylabel,
            title = "",
        )
        hspan!(
            fig[i],
            [target_lb[spec.index], target_ub[spec.index]];
            color = PLOT_COLORS[:target],
            alpha = 0.15,
            label = "target",
        )
    end

    savefig(fig, output_path)
    display(fig)
    return fig
end
