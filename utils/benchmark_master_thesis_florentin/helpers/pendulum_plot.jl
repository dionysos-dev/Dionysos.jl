if !isdefined(@__MODULE__, :_prepare_periodic_state_plot_data)
    include(joinpath(@__DIR__, "pendulum_plot_common.jl"))
end

function _pendulum_plot_filename(basename::AbstractString, stem::AbstractString)
    return _pendulum_artifact_filename(basename, stem)
end

function _pendulum_animation_filename(basename::AbstractString, stem::AbstractString, ext::AbstractString)
    return _pendulum_artifact_filename(basename, stem, ext)
end

function _prepare_pendulum_plot_data(
    problem,
    candidate;
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = SVector(1),
    periodic_periods = SVector(2pi),
    periodic_start = SVector(-pi),
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

function _pendulum_target_angle(problem)
    target_center = UT.get_center(problem.target_set)
    return target_center[1]
end

function _build_pendulum_pose_plot(
    xs,
    k::Int,
    target_angle::Real;
    title::AbstractString,
)
    θ = xs[k][1]
    bob_x = sin(θ)
    bob_y = -cos(θ)

    target_x = sin(target_angle)
    target_y = -cos(target_angle)

    fig = plot(
        xlim = (-1.25, 1.25),
        ylim = (-1.25, 0.45),
        aspect_ratio = :equal,
        legend = false,
        grid = false,
        framestyle = :box,
        title = title,
        xlabel = "x",
        ylabel = "y",
    )
    scatter!(fig, [0.0], [0.0]; color = :black, ms = 5)
    plot!(fig, [0.0, target_x], [0.0, target_y]; color = :red, lw = 2, ls = :dash)
    plot!(fig, [0.0, bob_x], [0.0, bob_y]; color = :royalblue4, lw = 4)
    scatter!(fig, [bob_x], [bob_y]; color = :royalblue4, ms = 8)
    return fig
end

function _build_pendulum_phase_plot(data, k::Int; title::AbstractString)
    return _build_phase_prefix_plot(data, k; dims = (1, 2), xlabel = L"\theta\,[\mathrm{rad}]", ylabel = L"\dot{\theta}\,[\mathrm{rad/s}]")
end

function _build_pendulum_control_plot(data, k::Int, problem; title::AbstractString)
    return _build_control_time_plot(data, problem; title = title, k = k)
end

function save_pendulum_animation!(
    output_dir::AbstractString,
    problem,
    candidate;
    basename::AbstractString = "",
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = SVector(1),
    periodic_periods = SVector(2pi),
    periodic_start = SVector(-pi),
    fps::Int = 15,
    every::Int = 1,
    save_gif::Bool = true,
    save_mp4::Bool = true,
    title::AbstractString = "Pendulum rollout",
)
    mkpath(output_dir)

    data = _prepare_pendulum_plot_data(
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

    target_angle = _pendulum_target_angle(problem)
    gif_path =
        save_gif ? joinpath(output_dir, _pendulum_animation_filename(basename, "rollout", "gif")) : nothing
    mp4_path =
        save_mp4 ? joinpath(output_dir, _pendulum_animation_filename(basename, "rollout", "mp4")) : nothing

    anim = @animate for k in 1:every:length(data.xs)
        fig_pose = _build_pendulum_pose_plot(data.xs, k, target_angle; title = title)
        fig_phase = _build_pendulum_phase_plot(data, k; title = "Phase space")
        fig_control = _build_pendulum_control_plot(data, max(k - 1, 1), problem; title = "Control")
        plot(fig_pose, fig_phase, fig_control; layout = @layout([a{0.42w} [b; c]]), size = (1200, 650))
    end

    if gif_path !== nothing
        gif(anim, gif_path; fps = fps)
    end

    if mp4_path !== nothing
        try
            mp4(anim, mp4_path; fps = fps)
        catch err
            @warn "Could not save pendulum MP4 animation." path = mp4_path exception = err
            mp4_path = nothing
        end
    end

    return (; gif_path, mp4_path)
end

function save_pendulum_plots!(
    output_dir::AbstractString,
    problem,
    candidate;
    basename::AbstractString = "",
    cert_result = nothing,
    show_ellipsoids::Bool = true,
    unwrap_angles::Bool = false,
    wrap_angles::Bool = true,
    periodic_dims = SVector(1),
    periodic_periods = SVector(2pi),
    periodic_start = SVector(-pi),
    phase_title::AbstractString = "Pendulum Phase Space",
    state_title::AbstractString = "Pendulum State vs Time",
    control_title::AbstractString = "Pendulum Control vs Time",
)
    mkpath(output_dir)

    data = _prepare_pendulum_plot_data(
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

    trajectory_label = occursin("mppi", lowercase(basename)) ? L"\mathrm{MPPI\ trajectory}" : L"\mathrm{nominal\ trajectory}"
    fig_phase = _build_state_space_projection(
        data;
        dims = (1, 2),
        xlabel = L"\theta\,[\mathrm{rad}]",
        ylabel = L"\dot{\theta}\,[\mathrm{rad/s}]",
        trajectory_label = trajectory_label,
        legend = :topright,
        size = (640, 460),
    )
    phase_path = joinpath(output_dir, _pendulum_plot_filename(basename, "phase_space"))
    savefig(fig_phase, phase_path)
    display(fig_phase)

    state_path = joinpath(output_dir, _pendulum_plot_filename(basename, "state_time"))
    _save_state_time_plot!(
        state_path,
        data,
        problem,
        (
            (; index = 1, label = L"\theta", ylabel = L"\theta\,[\mathrm{rad}]", color = :steelblue4),
            (; index = 2, label = L"\dot{\theta}", ylabel = L"\dot{\theta}\,[\mathrm{rad/s}]", color = :seagreen4),
        );
        title = "",
        layout = (2, 1),
        size = (900, 700),
    )

    fig_control = _build_control_time_plot(data, problem; title = "")
    control_path = joinpath(output_dir, _pendulum_plot_filename(basename, "control_time"))
    savefig(fig_control, control_path)
    display(fig_control)

    return (;
        phase_path,
        state_path,
        control_path,
    )
end
