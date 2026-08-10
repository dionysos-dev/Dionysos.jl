module CartPendulum

import LazySets
using StaticArrays
using MathematicalSystems
using Dionysos
using Plots

const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem

Base.@kwdef struct Params{T}
    M::T = 1.0      # cart mass
    m::T = 0.2      # pendulum mass
    l::T = 1.0      # distance to pendulum center of mass
    J::T = 0.0      # pendulum inertia around center of mass
    c::T = 0.1      # cart viscous friction
    γ::T = 0.05     # pendulum viscous friction
    g::T = 9.81
end

function dynamic(params::Params = Params())
    return (x, u) -> begin
        p = x[1]
        θ = x[2]
        v = x[3]
        ω = x[4]
        F = u[1]

        M = params.M
        m = params.m
        l = params.l
        J = params.J
        c = params.c
        γ = params.γ
        g = params.g

        Mt = M + m
        Jt = J + m*l^2

        sθ = sin(θ)
        cθ = cos(θ)

        A = @SMatrix [
            Mt -m*l*cθ
            -m*l*cθ Jt
        ]

        b = SVector(F - c*v - m*l*sθ*ω^2, m*g*l*sθ - γ*ω)

        acc = A \ b

        return SVector{4}(v, ω, acc[1], acc[2])
    end
end

"""
    jacobian_bound(params) -> u -> SMatrix{4,4}

Entrywise bound on `∂f/∂x`, as the growth-bound ODE `ṙ = L·r` requires: `L[i,j] ≥ |J[i,j]|`
off the diagonal and `L[i,i] ≥ J[i,i]` on it — the diagonal is signed, so a stabilising term
may legitimately be negative.

Obtained by bounding each symbolic Jacobian entry with interval arithmetic over the state set
of [`system`](@ref) and padding by 2%, the coupled cart–pendulum accelerations having no tidy
closed form. Passing no `jacobian_bound` makes Dionysos derive one the same way at solve time
(`Dionysos.System.compute_jacobian_bound`, needs `Symbolics` loaded).
"""
function jacobian_bound(params::Params = Params())
    return u -> SMatrix{4, 4}(
        # column 1 (cart position — the dynamics do not depend on it)
        0.0,
        0.0,
        0.0,
        0.0,
        # column 2 (pole angle)
        0.0,
        0.0,
        12.504,
        30.377,
        # column 3 (cart velocity)
        1.02,
        0.0,
        -0.081,
        0.103,
        # column 4 (pole angular velocity)
        0.0,
        1.02,
        2.092,
        1.735,
    )
end

function system(;
    params::Params = Params(),
    _X_ = LazySets.Hyperrectangle(;
        low = SVector(-5.0, -pi, -5.0, -5.0),
        high = SVector(5.0, pi, 5.0, 5.0),
    ),
    _U_ = LazySets.Hyperrectangle(; low = SVector(-10.0), high = SVector(10.0)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        LazySets.dim(_X_),
        LazySets.dim(_U_),
        _X_,
        _U_,
    )
end

function optimal_control_problem(;
    params::Params = Params(),
    objective = "swing_up",
    state_cost = nothing,
    transition_cost = nothing,
)
    if objective == "swing_up"
        _X_ = LazySets.Hyperrectangle(;
            low = SVector(-5.0, -pi, -5.0, -5.0),
            high = SVector(5.0, pi, 5.0, 5.0),
        )

        _U_ = LazySets.Hyperrectangle(; low = SVector(-10.0), high = SVector(10.0))

        _I_ = LazySets.Hyperrectangle(;
            low = SVector(-0.2, -5.0*pi/180.0, -0.2, -0.2),
            high = SVector(0.2, 5.0*pi/180.0, 0.2, 0.2),
        )

        _T_ = LazySets.Hyperrectangle(;
            low = SVector(-0.5, pi - 10.0*pi/180.0, -0.5, -0.5),
            high = SVector(0.5, pi + 10.0*pi/180.0, 0.5, 0.5),
        )
    elseif objective == "stabilize_up"
        _X_ = LazySets.Hyperrectangle(;
            low = SVector(-2.0, pi - 30.0*pi/180.0, -3.0, -3.0),
            high = SVector(2.0, pi + 30.0*pi/180.0, 3.0, 3.0),
        )

        _U_ = LazySets.Hyperrectangle(; low = SVector(-5.0), high = SVector(5.0))

        _I_ = LazySets.Hyperrectangle(;
            low = SVector(-0.5, pi - 15.0*pi/180.0, -0.5, -0.5),
            high = SVector(0.5, pi + 15.0*pi/180.0, 0.5, 0.5),
        )

        _T_ = LazySets.Hyperrectangle(;
            low = SVector(-0.2, pi - 5.0*pi/180.0, -0.2, -0.2),
            high = SVector(0.2, pi + 5.0*pi/180.0, 0.2, 0.2),
        )
    else
        error("Unknown objective: $objective")
    end

    sys = system(; params = params, _X_ = _X_, _U_ = _U_)

    return PR.OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost)
end

function system_plot!(;
    params::Params = Params(),
    cart_width = 0.6,
    cart_height = 0.3,
    pendulum_radius = 0.08,
    xlims = (-5.5, 5.5),
    ylims = (-1.4, 1.4),
    show_input = true,
)
    return function (fig, x, u)
        p = Float64(x[1])
        θ = Float64(x[2])
        F = Float64(u[1])
        l = Float64(params.l)

        cart_x = [
            p - cart_width / 2,
            p + cart_width / 2,
            p + cart_width / 2,
            p - cart_width / 2,
            p - cart_width / 2,
        ]

        cart_y = [
            -cart_height / 2,
            -cart_height / 2,
            cart_height / 2,
            cart_height / 2,
            -cart_height / 2,
        ]

        plot!(
            fig,
            cart_x,
            cart_y;
            linewidth = 2,
            fill = (true, 0.25),
            color = :blue,
            label = "",
        )

        pivot = SVector(p, cart_height / 2)
        bob = pivot + l * SVector(sin(θ), -cos(θ))

        plot!(
            fig,
            [pivot[1], bob[1]],
            [pivot[2], bob[2]];
            linewidth = 4,
            color = :black,
            label = "",
        )

        scatter!(
            fig,
            [pivot[1], bob[1]],
            [pivot[2], bob[2]];
            markersize = [5, 10],
            color = [:black, :red],
            label = "",
        )

        # ground rail
        plot!(
            fig,
            [xlims[1], xlims[2]],
            [-cart_height / 2, -cart_height / 2];
            linewidth = 2,
            color = :black,
            label = "",
        )

        if show_input
            annotate!(
                fig,
                xlims[1] + 0.05 * (xlims[2] - xlims[1]),
                ylims[1] + 0.08 * (ylims[2] - ylims[1]),
                text("F = $(round(F; digits = 2))", 10),
            )
        end

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end
