module DoublePendulum

import LazySets
using StaticArrays
using MathematicalSystems
using Plots
using Dionysos

const UT = Dionysos.Utils
const ST = Dionysos.System
const PB = Dionysos.Problem

Base.@kwdef struct Params{T}
    l1::T = 1.0
    l2::T = 1.0
    m1::T = 1.0
    m2::T = 1.0
    g::T = 9.81
end

function dynamic(params::Params = Params())
    return (x, u) -> begin
        M = params.m1 + params.m2
        Δθ = x[1] - x[2]
        α = params.m1 + params.m2 * sin(Δθ)^2

        return SVector{4}(
            x[3],
            x[4],
            (
                -sin(Δθ) * (
                    params.m2 * params.l1 * x[3]^2 * cos(Δθ) +
                    params.m2 * params.l2 * x[4]^2
                ) - params.g * (M * sin(x[1]) - params.m2 * sin(x[2]) * cos(Δθ))
            ) / (params.l1 * α) + u[1],
            (
                sin(Δθ) *
                (M * params.l1 * x[3]^2 + params.m2 * params.l2 * x[4]^2 * cos(Δθ)) +
                params.g * M * (sin(x[1]) * cos(Δθ) - sin(x[2]))
            ) / (params.l2 * α),
        )
    end
end

"""
    jacobian_bound(params) -> u -> SMatrix{4,4}

Entrywise bound on `∂f/∂x`, as the growth-bound ODE `ṙ = L·r` requires: `L[i,j] ≥ |J[i,j]|`
off the diagonal and `L[i,i] ≥ J[i,i]` on it.

The coupled accelerations have no tidy closed form, so this was obtained by bounding each
symbolic Jacobian entry with interval arithmetic over the state set of [`system`](@ref) and
padding by 2%. That domain contains the ones used by [`safety_problem`](@ref) and
[`optimal_control_problem`](@ref), and a bound valid on a set stays valid on its subsets, so
the same matrix serves all of them. It does not depend on `u`, which enters `f₃` additively.

Passing no `jacobian_bound` at all makes Dionysos derive one the same way at solve time
(`Dionysos.System.compute_jacobian_bound`, needs `Symbolics` loaded); this constant is the
cheaper, pre-computed equivalent.
"""
function jacobian_bound(params::Params = Params())
    return u -> SMatrix{4, 4}(
        # column 1 (θ₁)
        0.0,
        0.0,
        268.556,
        375.075,
        # column 2 (θ₂)
        0.0,
        0.0,
        258.55,
        375.075,
        # column 3 (ω₁)
        1.02,
        0.0,
        10.2,
        20.4,
        # column 4 (ω₂)
        0.0,
        1.02,
        10.2,
        10.2,
    )
end

function system(;
    params::Params = Params(),
    _X_ = UT.box(SVector(-π, -π, -5.0, -5.0), SVector(π, π, 5.0, 5.0)),
    _U_ = UT.box(SVector(-11.0), SVector(11.0)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        LazySets.dim(_X_),
        LazySets.dim(_U_),
        _X_,
        _U_,
    )
end

function safety_problem(; params::Params = Params(), objective = "safety_down")
    if objective == "safety_down"
        _X_ = UT.box(
            SVector(-π / 4.0, -π / 4.0, -5.0, -5.0),
            SVector(π / 4.0, π / 4.0, 5.0, 5.0),
        )

        _U_ = UT.box(SVector(-2.5), SVector(2.5))

        _I_ = UT.box(
            SVector(-3.0π / 180.0, -3.0π / 180.0, -0.5, -0.5),
            SVector(3.0π / 180.0, 3.0π / 180.0, 0.5, 0.5),
        )

        _S_ = UT.box(
            SVector(-35.0π / 180.0, -35.0π / 180.0, -1.0, -1.0),
            SVector(35.0π / 180.0, 35.0π / 180.0, 1.0, 1.0),
        )
    else
        error("Unknown objective: $objective")
    end

    sys = system(; params = params, _X_ = _X_, _U_ = _U_)
    return PB.SafetyProblem(sys, _I_, _S_)
end

function optimal_control_problem(;
    params::Params = Params(),
    objective = "swing_up",
    state_cost = nothing,
    transition_cost = nothing,
)
    if objective == "swing_up"
        _X_ = UT.box(SVector(-π / 2.0, -π, -5.0, -5.0), SVector(π / 2.0, π, 5.0, 5.0))

        _U_ = UT.set_minus(
            UT.box(SVector(-5.5), SVector(5.5)),
            UT.box(SVector(-0.5), SVector(0.5)),
        )

        _I_ = UT.box(
            SVector(-3.0π / 180.0, -3.0π / 180.0, -0.5, -0.5),
            SVector(3.0π / 180.0, 3.0π / 180.0, 0.5, 0.5),
        )

        _T_ = UT.box(
            SVector(-π / 2.0, π - 50.0π / 180.0, -4.5, -4.5),
            SVector(π / 2.0, π + 50.0π / 180.0, 4.5, 4.5),
        )
    else
        error("Unknown objective: $objective")
    end

    sys = system(; params = params, _X_ = _X_, _U_ = _U_)

    return PB.OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost)
end

function system_plot!(;
    params::Params = Params(),
    xlims = (-2.2, 2.2),
    ylims = (-2.2, 2.2),
    show_input = true,
    show_joints = true,
)
    return function (fig, x, u)
        θ1 = Float64(x[1])
        θ2 = Float64(x[2])
        τ = Float64(u[1])

        l1 = Float64(params.l1)
        l2 = Float64(params.l2)

        p0 = SVector(0.0, 0.0)
        p1 = p0 + l1 * SVector(sin(θ1), -cos(θ1))
        p2 = p1 + l2 * SVector(sin(θ2), -cos(θ2))

        plot!(
            fig,
            [p0[1], p1[1], p2[1]],
            [p0[2], p1[2], p2[2]];
            linewidth = 4,
            color = :black,
            label = "",
        )

        if show_joints
            scatter!(
                fig,
                [p0[1], p1[1], p2[1]],
                [p0[2], p1[2], p2[2]];
                markersize = [5, 8, 10],
                color = [:black, :blue, :red],
                label = "",
            )
        end

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
