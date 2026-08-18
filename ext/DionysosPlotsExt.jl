module DionysosPlotsExt

import Dionysos
const DI = Dionysos
const ST = DI.System

using Plots

function Dionysos.animate_trajectory_dashboard(
    system_plot!::Function,
    traj::ST.Trajectory;
    xdims = (1, 2),
    udims = (1,),
    Δt = 1.0, # physical time represented by each frame
    fps = 20,
    frame_step::Integer = 1, # render every `frame_step`-th frame; speeds up long trajectories
    filename::Union{Nothing, String} = nothing,
    title = "System evolution",
    xlabel_state = nothing,
    ylabel_state = nothing,
    xlabel_input = nothing,
    ylabel_input = nothing,
    xlims_system = nothing,
    ylims_system = nothing,
    xlims_state = nothing,
    ylims_state = nothing,
    xlims_input = nothing,
    ylims_input = nothing,
    show_full_state_traj = true,
    show_full_input_traj = true,
    # Optional `(p_state, x) -> p_state` background for the state panel, drawn
    # each frame before the trajectory — e.g. the slice of a carved state
    # region at the current state.
    state_background! = nothing,
    ylabel_mode = "mode",
)
    # The trajectory is self-describing: `states` is always present, the `inputs`
    # channel drives the input panel, and the `modes` channel is `nothing` for
    # non-hybrid systems (no mode panel) or the per-step mode otherwise (adds a
    # mode-vs-time panel, and `system_plot!` receives the mode as a 4th argument).
    ST.inputs(traj) === nothing &&
        error("animate_trajectory_dashboard needs a trajectory with inputs")
    xs = collect(ST.states(traj))
    us = collect(ST.inputs(traj))
    modevals = ST.modes(traj) === nothing ? nothing : collect(ST.modes(traj))

    N = length(xs)
    Nu = length(us)

    N == 0 && error("trajectory has no states")
    Nu == 0 && error("trajectory has no inputs")

    length(xdims) in (1, 2) || error("xdims must contain one or two dimensions")
    length(udims) in (1, 2) || error("udims must contain one or two dimensions")

    times_x = Δt .* collect(0:(N - 1))
    times_u = Δt .* collect(0:(Nu - 1))

    xvals = [[Float64(x[d]) for x in xs] for d in xdims]
    uvals = [[Float64(u[d]) for u in us] for d in udims]

    xlabel_state === nothing &&
        (xlabel_state = length(xdims) == 1 ? "time" : "x$(xdims[1])")
    ylabel_state === nothing &&
        (ylabel_state = length(xdims) == 1 ? "x$(xdims[1])" : "x$(xdims[2])")

    xlabel_input === nothing &&
        (xlabel_input = length(udims) == 1 ? "time" : "u$(udims[1])")
    ylabel_input === nothing &&
        (ylabel_input = length(udims) == 1 ? "u$(udims[1])" : "u$(udims[2])")

    # Render only every `frame_step`-th frame (always keeping the last) so long
    # trajectories animate fast. This is what actually cuts render time — `fps` only
    # sets playback speed, so all frames would otherwise still be drawn.
    frame_indices = frame_step <= 1 ? collect(1:N) : unique(vcat(1:frame_step:N, N))
    anim = @animate for k in frame_indices
        ku = min(k, Nu)

        # ----------------------------------------------------
        # Left: physical/system evolution
        # ----------------------------------------------------
        p_sys = plot(title = title, aspect_ratio = :equal, legend = false)

        xlims_system !== nothing && xlims!(p_sys, xlims_system...)
        ylims_system !== nothing && ylims!(p_sys, ylims_system...)

        if modevals === nothing
            system_plot!(p_sys, xs[k], us[ku])
        else
            system_plot!(p_sys, xs[k], us[ku], modevals[k])
        end

        # ----------------------------------------------------
        # Top-right: state evolution/projection
        # ----------------------------------------------------
        p_state = plot(
            xlabel = xlabel_state,
            ylabel = ylabel_state,
            title = length(xdims) == 1 ? "State evolution" : "State space",
            legend = :best,
        )

        xlims_state !== nothing && xlims!(p_state, xlims_state...)
        ylims_state !== nothing && ylims!(p_state, ylims_state...)

        state_background! !== nothing && state_background!(p_state, xs[k])

        if length(xdims) == 1
            if show_full_state_traj
                plot!(
                    p_state,
                    times_x,
                    xvals[1];
                    linestyle = :dot,
                    linewidth = 1,
                    label = "full state",
                )
            end

            plot!(p_state, times_x[1:k], xvals[1][1:k]; linewidth = 2, label = "past state")

            scatter!(
                p_state,
                [times_x[k]],
                [xvals[1][k]];
                markersize = 5,
                label = "current state",
            )

            vline!(p_state, [times_x[k]]; linestyle = :dash, label = "")
        else
            if show_full_state_traj
                plot!(
                    p_state,
                    xvals[1],
                    xvals[2];
                    linestyle = :dot,
                    linewidth = 1,
                    label = "full trajectory",
                )
            end

            plot!(
                p_state,
                xvals[1][1:k],
                xvals[2][1:k];
                linewidth = 2,
                label = "past trajectory",
            )

            scatter!(
                p_state,
                [xvals[1][k]],
                [xvals[2][k]];
                markersize = 5,
                label = "current state",
            )
        end

        # ----------------------------------------------------
        # Bottom-right: input evolution/projection
        # ----------------------------------------------------
        p_input = plot(
            xlabel = xlabel_input,
            ylabel = ylabel_input,
            title = length(udims) == 1 ? "Input evolution" : "Input space",
            legend = :best,
        )

        xlims_input !== nothing && xlims!(p_input, xlims_input...)
        ylims_input !== nothing && ylims!(p_input, ylims_input...)

        if length(udims) == 1
            if show_full_input_traj
                plot!(
                    p_input,
                    times_u,
                    uvals[1];
                    linestyle = :dot,
                    linewidth = 1,
                    label = "full input",
                )
            end

            plot!(
                p_input,
                times_u[1:ku],
                uvals[1][1:ku];
                linewidth = 2,
                label = "past input",
            )

            scatter!(
                p_input,
                [times_u[ku]],
                [uvals[1][ku]];
                markersize = 5,
                label = "current input",
            )

            vline!(p_input, [times_x[k]]; linestyle = :dash, label = "")
        else
            if show_full_input_traj
                plot!(
                    p_input,
                    uvals[1],
                    uvals[2];
                    linestyle = :dot,
                    linewidth = 1,
                    label = "full input trajectory",
                )
            end

            plot!(
                p_input,
                uvals[1][1:ku],
                uvals[2][1:ku];
                linewidth = 2,
                label = "past input trajectory",
            )

            scatter!(
                p_input,
                [uvals[1][ku]],
                [uvals[2][ku]];
                markersize = 5,
                label = "current input",
            )
        end

        # ----------------------------------------------------
        # Extra right panel: mode evolution (hybrid systems).
        # When present it leads the right column: mode, then state, then input.
        # ----------------------------------------------------
        if modevals === nothing
            panels = Any[p_sys, p_state, p_input]
        else
            p_mode = plot(;
                xlabel = "time",
                ylabel = ylabel_mode,
                title = "Mode evolution",
                legend = false,
                yticks = sort(unique(modevals)),
            )
            plot!(
                p_mode,
                times_x,
                modevals;
                seriestype = :steppost,
                linewidth = 1,
                linestyle = :dot,
            )
            plot!(
                p_mode,
                times_x[1:k],
                modevals[1:k];
                seriestype = :steppost,
                linewidth = 2,
            )
            scatter!(p_mode, [times_x[k]], [modevals[k]]; markersize = 5)
            vline!(p_mode, [times_x[k]]; linestyle = :dash)
            panels = Any[p_sys, p_mode, p_state, p_input]
        end

        layout =
            modevals === nothing ? @layout([a{0.48w} [b; c]]) :
            @layout([a{0.48w} [b; c; d]])
        plot(
            panels...;
            layout = layout,
            size = (1400, 800),
            plot_title = "t = $(round(times_x[k]; digits = 2)) s   frame $k / $N",
        )
    end

    if filename === nothing
        return anim
    end
    if endswith(lowercase(filename), ".gif")
        gif(anim, filename; fps = fps)
    elseif endswith(lowercase(filename), ".mp4")
        mp4(anim, filename; fps = fps)
    else
        error("filename must end with .gif or .mp4")
    end

    return filename
end

end
