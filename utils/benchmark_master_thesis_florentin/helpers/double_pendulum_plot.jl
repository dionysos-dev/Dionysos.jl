if !isdefined(@__MODULE__, :_prepare_periodic_state_plot_data)
    include(joinpath(@__DIR__, "pendulum_plot_common.jl"))
end

function _double_pendulum_plot_filename(basename::AbstractString, stem::AbstractString)
    return _pendulum_artifact_filename(basename, stem)
end

function _double_pendulum_animation_filename(
    basename::AbstractString,
    stem::AbstractString,
    ext::AbstractString,
)
    return _pendulum_artifact_filename(basename, stem, ext)
end

function _prepare_double_pendulum_plot_data(
    problem,
    candidate;
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = SVector(1, 2),
    periodic_periods = SVector(2pi, 2pi),
    periodic_start = SVector(-pi, -pi),
)
    return _prepare_periodic_state_plot_data(
        problem,
        candidate;
        cert_result = cert_result,
        show_ellipsoids = show_ellipsoids,
        unwrap_angles = unwrap_angles,
        wrap_angles = wrap_angles,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )
end

function _double_pendulum_target_angles(problem)
    target_center = UT.get_center(problem.target_set)
    return target_center[1], target_center[2]
end

function _build_double_pendulum_pose_plot(
    xs,
    k::Int,
    target_angles;
    l1::Float64,
    l2::Float64,
    title::AbstractString,
)
    θ1 = xs[k][1]
    θ2 = xs[k][2]
    θ1t, θ2t = target_angles

    x1 = l1 * sin(θ1)
    y1 = -l1 * cos(θ1)
    x2 = x1 + l2 * sin(θ2)
    y2 = y1 - l2 * cos(θ2)

    xt1 = l1 * sin(θ1t)
    yt1 = -l1 * cos(θ1t)
    xt2 = xt1 + l2 * sin(θ2t)
    yt2 = yt1 - l2 * cos(θ2t)

    reach = l1 + l2 + 0.15
    fig = plot(
        xlim = (-reach, reach),
        ylim = (-reach, reach),
        aspect_ratio = :equal,
        legend = false,
        grid = false,
        framestyle = :box,
        title = title,
        xlabel = "x",
        ylabel = "y",
    )
    scatter!(fig, [0.0], [0.0]; color = :black, ms = 5)
    plot!(fig, [0.0, xt1, xt2], [0.0, yt1, yt2]; color = :red, lw = 2, ls = :dash)
    plot!(fig, [0.0, x1, x2], [0.0, y1, y2]; color = :royalblue4, lw = 4)
    scatter!(fig, [x1, x2], [y1, y2]; color = :royalblue4, ms = 7)
    return fig
end

function _build_double_pendulum_angles_plot(data, k::Int; title::AbstractString)
    return _build_phase_prefix_plot(data, k; dims = (1, 2), xlabel = "θ₁ [rad]", ylabel = "θ₂ [rad]")
end

function _build_double_pendulum_control_plot(data, k::Int, problem; title::AbstractString)
    return _build_control_time_plot(data, problem; title = title, k = k)
end

function save_double_pendulum_animation!(
    output_dir::AbstractString,
    problem,
    candidate;
    basename::AbstractString = "",
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = SVector(1, 2),
    periodic_periods = SVector(2pi, 2pi),
    periodic_start = SVector(-pi, -pi),
    l1::Float64 = 1.0,
    l2::Float64 = 1.0,
    fps::Int = 15,
    every::Int = 1,
    save_gif::Bool = true,
    save_mp4::Bool = true,
    title::AbstractString = "Double pendulum rollout",
)
    mkpath(output_dir)

    data = _prepare_double_pendulum_plot_data(
        problem,
        candidate;
        cert_result = cert_result,
        show_ellipsoids = show_ellipsoids,
        unwrap_angles = unwrap_angles,
        wrap_angles = wrap_angles,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )

    target_angles = _double_pendulum_target_angles(problem)
    gif_path = save_gif ?
        joinpath(output_dir, _double_pendulum_animation_filename(basename, "rollout", "gif")) :
        nothing
    mp4_path = save_mp4 ?
        joinpath(output_dir, _double_pendulum_animation_filename(basename, "rollout", "mp4")) :
        nothing

    anim = @animate for k in 1:every:length(data.xs)
        fig_pose = _build_double_pendulum_pose_plot(
            data.xs,
            k,
            target_angles;
            l1 = l1,
            l2 = l2,
            title = title,
        )
        fig_angles = _build_double_pendulum_angles_plot(data, k; title = "(theta1, theta2)")
        fig_control = _build_double_pendulum_control_plot(data, max(k - 1, 1), problem; title = "Control")
        plot(fig_pose, fig_angles, fig_control; layout = @layout([a{0.42w} [b; c]]), size = (1200, 650))
    end

    if gif_path !== nothing
        gif(anim, gif_path; fps = fps)
    end

    if mp4_path !== nothing
        try
            mp4(anim, mp4_path; fps = fps)
        catch err
            @warn "Could not save double pendulum MP4 animation." path = mp4_path exception = err
            mp4_path = nothing
        end
    end

    return (; gif_path, mp4_path)
end

function save_double_pendulum_plots!(
    output_dir::AbstractString,
    problem,
    candidate;
    basename::AbstractString = "",
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = SVector(1, 2),
    periodic_periods = SVector(2pi, 2pi),
    periodic_start = SVector(-pi, -pi),
    angles_title::AbstractString = "Double pendulum angles",
    velocities_title::AbstractString = "Double pendulum velocities",
    state_title::AbstractString = "Double pendulum states",
    control_title::AbstractString = "Double pendulum control",
    phase_title::AbstractString = "Double pendulum phase portraits",
)
    mkpath(output_dir)

    data = _prepare_double_pendulum_plot_data(
        problem,
        candidate;
        cert_result = cert_result,
        show_ellipsoids = show_ellipsoids,
        unwrap_angles = unwrap_angles,
        wrap_angles = wrap_angles,
        periodic_dims = periodic_dims,
        periodic_periods = periodic_periods,
        periodic_start = periodic_start,
    )

    trajectory_label = occursin("mppi", lowercase(basename)) ? "MPPI trajectory" : "nominal trajectory"
    fig_angles = _build_state_space_projection(
        data;
        dims = (1, 2),
        xlabel = "θ₁ [rad]",
        ylabel = "θ₂ [rad]",
        trajectory_label = trajectory_label,
    )
    angles_path = joinpath(output_dir, _double_pendulum_plot_filename(basename, "angles_space"))
    savefig(fig_angles, angles_path)
    display(fig_angles)

    fig_velocities = _build_state_space_projection(
        data;
        dims = (3, 4),
        use_raw_sets = true,
        xlabel = "θ̇₁ [rad/s]",
        ylabel = "θ̇₂ [rad/s]",
        trajectory_label = trajectory_label,
    )
    velocities_path =
        joinpath(output_dir, _double_pendulum_plot_filename(basename, "velocity_space"))
    savefig(fig_velocities, velocities_path)
    display(fig_velocities)

    state_path = joinpath(output_dir, _double_pendulum_plot_filename(basename, "state_time"))
    _save_state_time_plot!(
        state_path,
        data,
        problem,
        (
            (; index = 1, label = "θ₁", ylabel = "θ₁ [rad]", color = :royalblue4),
            (; index = 2, label = "θ₂", ylabel = "θ₂ [rad]", color = :darkorange3),
            (; index = 3, label = "θ̇₁", ylabel = "θ̇₁ [rad/s]", color = :darkgreen),
            (; index = 4, label = "θ̇₂", ylabel = "θ̇₂ [rad/s]", color = :purple4),
        );
        title = "",
        layout = (2, 2),
        size = (1100, 800),
    )

    fig_phase = plot(; layout = (1, 2), thesis_plot_kwargs(; legend = :topright, size = (1100, 450))...)
    sets = _plot_data_sets(data)
    plot_state_space_basic!(
        fig_phase[1],
        sets.domain,
        sets.initial_set,
        sets.target_set,
        data.x_traj;
        dims = [1, 3],
        trajectory_label = trajectory_label,
    )
    _add_projected_ellipsoids!(fig_phase[1], data.ellipsoids, (1, 3))
    plot!(fig_phase[1]; title = "", xlabel = "θ₁ [rad]", ylabel = "θ̇₁ [rad/s]")
    plot_state_space_basic!(
        fig_phase[2],
        sets.domain,
        sets.initial_set,
        sets.target_set,
        data.x_traj;
        dims = [2, 4],
        trajectory_label = trajectory_label,
    )
    _add_projected_ellipsoids!(fig_phase[2], data.ellipsoids, (2, 4))
    plot!(fig_phase[2]; title = "", xlabel = "θ₂ [rad]", ylabel = "θ̇₂ [rad/s]")
    phase_path = joinpath(output_dir, _double_pendulum_plot_filename(basename, "phase_portraits"))
    savefig(fig_phase, phase_path)
    display(fig_phase)

    fig_control = _build_control_time_plot(data, problem; title = "")
    control_path = joinpath(output_dir, _double_pendulum_plot_filename(basename, "control_time"))
    savefig(fig_control, control_path)
    display(fig_control)

    return (;
        angles_path,
        velocities_path,
        state_path,
        phase_path,
        control_path,
    )
end
