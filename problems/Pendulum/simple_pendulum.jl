module SimplePendulum

import LazySets
using StaticArrays
using MathematicalSystems
using Plots
using Dionysos

const UT = Dionysos.Utils
const PR = Dionysos.Problem
const ST = Dionysos.System

Base.@kwdef struct Params{T}
    l::T = 1.0
    g::T = 9.81
end

function dynamic(params::Params = Params())
    return (x, u) -> SVector{2}(x[2], -(params.g / params.l) * sin(x[1]) + u[1])
end

# The growth bound drives `ṙ = L * r`, so `L` must dominate the Jacobian entrywise:
#
#     J = [0                  1]        |J| ≤ L = [0      1]
#         [-(g/l)cos(x₁)      0]                  [g/l    0]
#
# `SMatrix` fills column by column, so the second and third arguments are the *lower-left* and
# *upper-right* entries — transposing them here silently under-bounds the velocity radius and
# the abstraction stops covering the real dynamics.
function jacobian_bound(params::Params = Params())
    return u -> SMatrix{2, 2}(0.0, params.g / params.l, 1.0, 0.0)
end

function system(;
    params::Params = Params(),
    _X_ = UT.box(SVector(-π, -5.0), SVector(π, 5.0)),
    _U_ = UT.box(SVector(-6.0), SVector(6.0)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        LazySets.dim(_X_),
        LazySets.dim(_U_),
        _X_,
        _U_,
    )
end

function safety_problem(; params::Params = Params(), objective = "safety_up")
    if objective == "safety_down"
        _X_ = UT.box(SVector(-π, -1.5), SVector(π, 1.5))
        _U_ = UT.box(SVector(-4.0), SVector(4.0))
        _I_ = UT.box(SVector(-3.0π / 180.0, -0.5), SVector(3.0π / 180.0, 0.5))
        _S_ = UT.box(SVector(-15.0π / 180.0, -1.0), SVector(15.0π / 180.0, 1.0))
    elseif objective == "safety_up"
        _X_ = UT.box(SVector(-π, -1.5), SVector(π, 1.5))
        _U_ = UT.box(SVector(-4.0), SVector(4.0))
        _I_ = UT.box(SVector(π - 3.0π / 180.0, -0.5), SVector(π + 3.0π / 180.0, 0.5))
        _S_ = UT.box(SVector(π - 15.0π / 180.0, -1.0), SVector(π + 15.0π / 180.0, 1.0))
    else
        error("Unknown objective: $objective")
    end

    sys = system(; params = params, _X_ = _X_, _U_ = _U_)
    return PR.SafetyProblem(sys, _I_, _S_)
end

function optimal_control_problem(;
    params::Params = Params(),
    objective = "reachability_up_low_power",
    state_cost = nothing,
    transition_cost = nothing,
)
    _O_ = nothing

    if objective == "reachability_up_high_power"
        _X_ = UT.box(SVector(-π, -5.0), SVector(π, 5.0))
        _U_ = UT.set_minus(
            UT.box(SVector(-10.0), SVector(10.0)),
            UT.box(SVector(-0.5), SVector(0.5)),
        )
        _I_ = UT.box(SVector(-5.0π / 180.0, -0.2), SVector(5.0π / 180.0, 0.2))
        _T_ = UT.box(SVector(π - 15.0π / 180.0, -1.0), SVector(π + 15.0π / 180.0, 1.0))

    elseif objective == "reachability_up_medium_power"
        _X_ = UT.box(SVector(-π, -5.0), SVector(π, 5.0))
        _U_ = UT.set_minus(
            UT.box(SVector(-7.0), SVector(7.0)),
            UT.box(SVector(-0.5), SVector(0.5)),
        )
        _I_ = UT.box(SVector(-5.0π / 180.0, -0.2), SVector(5.0π / 180.0, 0.2))
        _T_ = UT.box(SVector(π - 15.0π / 180.0, -1.0), SVector(π + 15.0π / 180.0, 1.0))

    elseif objective == "reachability_up_medium_power_no_obstacle"
        _X_ = UT.box(SVector(-π, -7.0), SVector(π, 7.0))
        _U_ = UT.box(SVector(-4.5), SVector(4.5))
        _I_ = UT.box(SVector(-10.0π / 180.0, -0.5), SVector(10.0π / 180.0, 0.5))
        _T_ = UT.box(SVector(π - 15.0π / 180.0, -1.0), SVector(π + 15.0π / 180.0, 1.0))

    elseif objective == "reachability_up_low_power"
        _X_ = UT.box(SVector(-π, -7.0), SVector(π, 7.0))
        _U_ = UT.set_minus(
            UT.box(SVector(-2.5), SVector(2.5)),
            UT.box(SVector(-0.5), SVector(0.5)),
        )
        _I_ = UT.box(SVector(-5.0π / 180.0, -0.2), SVector(5.0π / 180.0, 0.2))
        _T_ = UT.box(SVector(π - 15.0π / 180.0, -1.0), SVector(π + 15.0π / 180.0, 1.0))
        _O_ = UT.box(SVector(-π + 16.0π / 180.0, -7.0), SVector(-π + 38.0π / 180.0, 7.0))

    else
        error("Unknown objective: $objective")
    end

    _X_ = _O_ === nothing ? _X_ : UT.set_minus(_X_, _O_)
    sys = system(; params = params, _X_ = _X_, _U_ = _U_)

    return PR.OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost)
end

function system_plot!(;
    params::Params = Params(),
    xlims = (-1.3params.l, 1.3params.l),
    ylims = (-1.3params.l, 1.3params.l),
    show_input = true,
    show_pivot = true,
)
    return function (fig, x, u)
        θ = Float64(x[1])
        τ = Float64(u[1])
        l = Float64(params.l)

        px = l * sin(θ)
        py = -l * cos(θ)

        if show_pivot
            scatter!(fig, [0.0], [0.0]; markersize = 5, color = :black, label = "")
        end

        plot!(fig, [0.0, px], [0.0, py]; linewidth = 4, color = :black, label = "")

        scatter!(fig, [px], [py]; markersize = 10, color = :blue, label = "")

        if show_input
            annotate!(
                fig,
                xlims[1] + 0.05 * (xlims[2] - xlims[1]),
                ylims[1] + 0.05 * (ylims[2] - ylims[1]),
                text("τ = $(round(τ; digits = 2))", 10),
            )
        end

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end
