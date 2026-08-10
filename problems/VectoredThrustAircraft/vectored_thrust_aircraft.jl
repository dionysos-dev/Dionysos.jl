module VectoredThrustAircraft
# The model is nonlinear, underactuated-looking in position/orientation, and it has a clean 2D physical animation.
import LazySets
using StaticArrays
using MathematicalSystems
using Plots

using Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem

Base.@kwdef struct Params{T}
    m::T = 1.0      # mass
    J::T = 1.0      # moment of inertia
    r::T = 1.0      # thrust offset
    g::T = 9.81     # gravity
    c::T = 0.1      # translational damping
end

# State:
# x = [px, py, θ, vx, vy, ω]
#
# Input:
# u = [u1, u2]
#
# Book convention:
# u1 = F1
# u2 = F2 - m*g
#
# Dynamics from Example 2.9:
# m ẍ = -m g sin θ - c ẋ + u1 cos θ - u2 sin θ
# m ÿ =  m g (cos θ - 1) - c ẏ + u1 sin θ + u2 cos θ
# J θ̈ = r u1
function dynamic(params::Params = Params())
    return (x, u) -> begin
        px = x[1] # horizontal position
        py = x[2] # vertical position
        θ = x[3] # attitude (pitch angle) Attitude: how the aircraft is rotated
        vx = x[4] # horizontal velocity
        vy = x[5] # vertical velocity
        ω = x[6] # angular velocity

        u1 = u[1] # thrust along body x-axis
        u2 = u[2] # thrust along body y-axis (in the simplified model, this is a lateral force)

        return SVector{6}(
            vx,
            vy,
            ω,
            -params.g * sin(θ) - (params.c / params.m) * vx + (u1 / params.m) * cos(θ) -
            (u2 / params.m) * sin(θ),
            params.g * (cos(θ) - 1.0) - (params.c / params.m) * vy +
            (u1 / params.m) * sin(θ) +
            (u2 / params.m) * cos(θ),
            (params.r / params.J) * u1,
        )
    end
end

function jacobian_bound(params::Params = Params())
    return u -> begin
        u1 = abs(u[1])
        u2 = abs(u[2])

        # Conservative bound on |∂f/∂x|.
        #
        # f4 depends on θ and vx:
        # |∂f4/∂θ| ≤ g + |u1|/m + |u2|/m
        # |∂f4/∂vx| = c/m
        #
        # f5 depends on θ and vy:
        # |∂f5/∂θ| ≤ g + |u1|/m + |u2|/m
        # |∂f5/∂vy| = c/m
        bθ = params.g + (u1 + u2) / params.m
        bv = params.c / params.m

        # `SMatrix` fills column by column, so each group of six below is a *column*:
        #
        #   row 1..3 : ẋ = vx, ẏ = vy, θ̇ = ω          → 1 in columns 4, 5, 6
        #   row 4    : bθ in column 3 (θ), -bv in column 4 (vx)
        #   row 5    : bθ in column 3 (θ), -bv in column 5 (vy)
        #
        # The velocity terms sit on the diagonal, where the bound is on the *signed* entry
        # rather than its magnitude — drag is stabilising, so -c/m is both correct and
        # tighter than +c/m.
        return SMatrix{6, 6}(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            bθ,
            bθ,
            0.0,
            1.0,
            0.0,
            0.0,
            -bv,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            -bv,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
        )
    end
end

function system(;
    params::Params = Params(),
    _X_ = LazySets.Hyperrectangle(;
        low = SVector(-5.0, -5.0, -pi / 2, -4.0, -4.0, -4.0),
        high = SVector(5.0, 5.0, pi / 2, 4.0, 4.0, 4.0),
    ),
    _U_ = LazySets.Hyperrectangle(; low = SVector(-5.0, -5.0), high = SVector(5.0, 5.0)),
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
    _X_ = LazySets.Hyperrectangle(;
        low = SVector(-5.0, -5.0, -pi / 2, -4.0, -4.0, -4.0),
        high = SVector(5.0, 5.0, pi / 2, 4.0, 4.0, 4.0),
    ),
    _U_ = LazySets.Hyperrectangle(; low = SVector(-5.0, -5.0), high = SVector(5.0, 5.0)),
    _I_ = LazySets.Hyperrectangle(;
        low = SVector(-4.5, -4.5, -0.1, -0.2, -0.2, -0.2),
        high = SVector(-4.0, -4.0, 0.1, 0.2, 0.2, 0.2),
    ),
    _T_ = LazySets.Hyperrectangle(;
        low = SVector(4.0, 4.0, -0.1, -0.3, -0.3, -0.3),
        high = SVector(4.5, 4.5, 0.1, 0.3, 0.3, 0.3),
    ),
    state_cost = nothing,
    transition_cost = nothing,
)
    sys = system(; params = params, _X_ = _X_, _U_ = _U_)

    return PR.OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost)
end

function system_plot!(;
    body_length = 0.8,
    body_width = 0.35,
    thrust_scale = 0.15,
    xlims = (-6.0, 6.0),
    ylims = (-6.0, 6.0),
    show_input = true,
)
    return function (fig, x, u)
        px = Float64(x[1])
        py = Float64(x[2])
        θ = Float64(x[3])

        u1 = Float64(u[1])
        u2 = Float64(u[2])

        R = @SMatrix [
            cos(θ) -sin(θ)
            sin(θ) cos(θ)
        ]

        c0 = SVector(px, py)

        body_local = [
            SVector(body_length / 2, body_width / 2),
            SVector(body_length / 2, -body_width / 2),
            SVector(-body_length / 2, -body_width / 2),
            SVector(-body_length / 2, body_width / 2),
            SVector(body_length / 2, body_width / 2),
        ]

        body = [c0 + R * q for q in body_local]

        plot!(
            fig,
            [p[1] for p in body],
            [p[2] for p in body];
            linewidth = 2,
            fill = (true, 0.25),
            color = :blue,
            label = "",
        )

        # Body heading
        nose = c0 + R * SVector(body_length / 2, 0.0)
        plot!(fig, [px, nose[1]], [py, nose[2]]; linewidth = 3, color = :blue, label = "")

        # Draw thrust vector in body coordinates.
        # u1 acts along body x-axis, u2 along body y-axis in the simplified model.
        thrust_body = SVector(u1, u2)
        thrust_world = R * thrust_body

        p_force = c0 - thrust_scale * thrust_world

        plot!(
            fig,
            [px, p_force[1]],
            [py, p_force[2]];
            linewidth = 3,
            color = :red,
            linestyle = :dash,
            label = "",
        )

        scatter!(fig, [px], [py]; markersize = 5, color = :black, label = "")

        if show_input
            annotate!(
                fig,
                xlims[1] + 0.05 * (xlims[2] - xlims[1]),
                ylims[1] + 0.05 * (ylims[2] - ylims[1]),
                text(
                    "u₁ = $(round(u1; digits = 2))\n" * "u₂ = $(round(u2; digits = 2))",
                    10,
                ),
            )
        end

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end
