module UnicycleRobot

using StaticArrays
using Plots

"""
    system_plot!(; xlims, ylims, arrow_length = 0.5, markersize = 7, color = :blue)

Return the `(fig, x, u) -> fig` drawing closure for the unicycle cart: a marker at the
position `(x₁, x₂)` with a stick pointing along the heading `x₃`. Hand it to
`Dionysos.animate_trajectory_dashboard` as the system view.
"""
function system_plot!(; xlims, ylims, arrow_length = 0.5, markersize = 7, color = :blue)
    return function (fig, x, u)
        θ = Float64(x[3])
        c = SVector(Float64(x[1]), Float64(x[2]))
        head = c + arrow_length * SVector(cos(θ), sin(θ))

        scatter!(fig, [c[1]], [c[2]]; markersize = markersize, color = color, label = "")
        plot!(fig, [c[1], head[1]], [c[2], head[2]]; lw = 3, color = color, label = "")

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end
