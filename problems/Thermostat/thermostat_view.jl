# Shared thermostat system view. Included inside each `Thermostat*` module so the
# continuous and the two hybrid variants render an identical thermometer in the
# `animate_trajectory_dashboard` system panel. Kept here (rather than copied per
# module) so the three views stay aligned.
#
# The including module provides `using Plots` in scope.

const _THERMO_TLIMS = (16.0, 26.0)

# Draw the thermostat view on `fig`: a mercury column filled to temperature `T`
# (orange when `heater_on`, blue otherwise) over the target band `[tlo, thi]`, plus a
# top annotation. `aspect_ratio = :auto` overrides the dashboard's square system panel
# (meant for 2-D views) so this 1-D thermometer fills the panel width.
function _draw_thermometer!(fig, T, heater_on, tlo, thi; tlims = _THERMO_TLIMS, label = "")
    col = heater_on ? :orangered : :steelblue

    # Target temperature band (spans the full width).
    plot!(
        fig,
        [0.0, 1.0, 1.0, 0.0, 0.0],
        [tlo, tlo, thi, thi, tlo];
        aspect_ratio = :auto,
        seriestype = :shape,
        fillcolor = :seagreen,
        fillalpha = 0.15,
        linecolor = :seagreen,
        linestyle = :dash,
        label = "",
    )

    # Mercury column filled up to the current temperature.
    plot!(
        fig,
        [0.35, 0.65, 0.65, 0.35, 0.35],
        [tlims[1], tlims[1], T, T, tlims[1]];
        seriestype = :shape,
        fillcolor = col,
        linecolor = col,
        label = "",
    )
    scatter!(fig, [0.5], [T]; color = col, markersize = 14, label = "")

    annotate!(
        fig,
        0.5,
        tlims[2] - 0.4,
        text((isempty(label) ? "" : "$label\n") * "T = $(round(T; digits = 2)) °C", 9),
    )

    xlims!(fig, 0.0, 1.0)
    ylims!(fig, tlims...)
    return fig
end
