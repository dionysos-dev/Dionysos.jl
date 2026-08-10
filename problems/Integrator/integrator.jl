module Integrator

import LazySets
using StaticArrays
using MathematicalSystems
using Plots

using Dionysos

const UT = Dionysos.Utils
const PR = Dionysos.Problem

# 2D integrator: ẋ = u
dynamic() = (x, u) -> SVector{2}(u[1], u[2])

# For ẋ = u, ∂f/∂x = 0, so a constant zero matrix bound is valid.
# Returned function signature matches Dionysos' growth-bound constructors: L(u).
function jacobian_bound()
    return (u) -> @SMatrix [
        0.0 0.0;
        0.0 0.0
    ]
end

function system(;
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0)),
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(),
        LazySets.dim(_X_),
        LazySets.dim(_U_),
        _X_,
        _U_,
    )
end

function optimal_control_problem(;
    _X_ = LazySets.Hyperrectangle(; low = SVector(-2.0, -2.0), high = SVector(2.0, 2.0)),
    _U_ = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0)),
    _I_ = LazySets.Hyperrectangle(; low = SVector(-1.6, -1.6), high = SVector(-1.4, -1.4)),
    _T_ = LazySets.Hyperrectangle(; low = SVector(-0.2, -0.2), high = SVector(0.2, 0.2)),
)
    sys = system(; _X_ = _X_, _U_ = _U_)
    return PR.OptimalControlProblem(sys, _I_, _T_, nothing, nothing)
end

function system_plot!(;
    xlims = (-2.2, 2.2),
    ylims = (-2.2, 2.2),
    marker_size = 8,
    show_input = true,
)
    return function (fig, x, u)
        px = Float64(x[1])
        py = Float64(x[2])
        ux = Float64(u[1])
        uy = Float64(u[2])

        scatter!(fig, [px], [py]; markersize = marker_size, color = :blue, label = "")

        plot!(
            fig,
            [px, px + ux],
            [py, py + uy];
            linewidth = 3,
            color = :red,
            linestyle = :dash,
            label = "",
        )

        if show_input
            annotate!(
                fig,
                xlims[1] + 0.05 * (xlims[2] - xlims[1]),
                ylims[1] + 0.05 * (ylims[2] - ylims[1]),
                text(
                    "u₁ = $(round(ux; digits = 2))\n" * "u₂ = $(round(uy; digits = 2))",
                    10,
                ),
            )
        end

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end # module
