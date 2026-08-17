module ArticulatedVehicle2Trailers

using StaticArrays
using MathematicalSystems
using Dionysos
import LazySets
using Plots

const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem

# Tractor towing TWO trailers (ported from PR #569's 5-D model, on the modern
# API and the sibling 4-D file's input convention u = [v, δ] — the steering
# nonlinearity tan δ lives inside the dynamics, not in the input coordinate).
#
# State x = [x1, x2, θ1, ϕ1, ϕ2]: tractor rear-axle position, tractor heading,
# hitch angle tractor→trailer 1, hitch angle trailer 1→trailer 2. Body
# headings: θ2 = θ1 + ϕ1, θ3 = θ2 + ϕ2.

Base.@kwdef struct Params{T}
    L1::T = 2.2      # tractor wheelbase
    L2::T = 2.8      # hitch 1 → trailer-1 axle
    L3::T = 2.2      # hitch 2 → trailer-2 axle
    Lc::T = 0.9      # off-axle hitch offset behind the tractor
    Lc2::T = 0.9     # off-axle hitch offset behind trailer 1
end

# ----------------------------
# Dynamics: ẋ = f(x, u), u = [v, δ]
# ----------------------------
function dynamic(p::Params = Params())
    return (x, u) -> begin
        v = u[1]
        tδ = tan(u[2])
        θ1 = x[3]
        ϕ1 = x[4]
        ϕ2 = x[5]

        return SVector{5}(
            v * cos(θ1),
            v * sin(θ1),
            (v / p.L1) * tδ,
            -(v / (p.L1 * p.L2)) * (p.L1 * sin(ϕ1) + (p.Lc * cos(ϕ1) + p.L2) * tδ),
            -(v / p.L3) * (cos(ϕ1) - (p.Lc / p.L1) * tδ * sin(ϕ1)) * sin(ϕ2) +
            (1 + (p.Lc2 / p.L3) * cos(ϕ2)) *
            ((v / p.L2) * sin(ϕ1) + ((v * p.Lc) / (p.L1 * p.L2)) * tδ * cos(ϕ1)),
        )
    end
end

# ----------------------------
# Jacobian A(x, u) = ∂f/∂x — only the θ1 and ϕ columns are nonzero.
# ----------------------------
function jacobian(p::Params = Params())
    return (x, u) -> begin
        v = u[1]
        tδ = tan(u[2])
        θ1 = x[3]
        ϕ1 = x[4]
        ϕ2 = x[5]

        A = cos(ϕ1) - (p.Lc / p.L1) * tδ * sin(ϕ1)
        B = (1 / p.L2) * sin(ϕ1) + (p.Lc / (p.L1 * p.L2)) * tδ * cos(ϕ1)
        Aϕ = -sin(ϕ1) - (p.Lc / p.L1) * tδ * cos(ϕ1)
        Bϕ = (1 / p.L2) * cos(ϕ1) - (p.Lc / (p.L1 * p.L2)) * tδ * sin(ϕ1)

        d4dϕ1 = -(v / (p.L1 * p.L2)) * (p.L1 * cos(ϕ1) - p.Lc * tδ * sin(ϕ1))
        d5dϕ1 = -(v / p.L3) * Aϕ * sin(ϕ2) + v * (1 + (p.Lc2 / p.L3) * cos(ϕ2)) * Bϕ
        d5dϕ2 = -(v / p.L3) * A * cos(ϕ2) - v * (p.Lc2 / p.L3) * B * sin(ϕ2)

        return @SMatrix [
            0.0 0.0 -v*sin(θ1) 0.0 0.0
            0.0 0.0 v*cos(θ1) 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 d4dϕ1 0.0
            0.0 0.0 0.0 d5dϕ1 d5dϕ2
        ]
    end
end

# ----------------------------
# Elementwise Jacobian bound |A(x, u)| ≤ Jb(u): harmonic-amplitude form.
# Every ϕ-entry is a linear combination of sin/cos with coefficients built
# from R = √(1 + ((Lc/L1)·tanδ)²): |A|, |Aϕ| ≤ R and |B|, |Bϕ| ≤ R/L2, so
#   |∂f4/∂ϕ1| ≤ (v/L2)·R
#   |∂f5/∂ϕ1| ≤ v·(R/L3 + (1 + |Lc2|/L3)·R/L2)
#   |∂f5/∂ϕ2| ≤ v·(R/L3 + (|Lc2|/L3)·R/L2)
# — strictly tighter than the triangle bounds (PR #569's one real algorithmic
# improvement in this file, verified against the exact Jacobian).
# ----------------------------
function jacobian_bound(p::Params = Params())
    return u -> begin
        v = abs(u[1])
        tδ = abs(tan(u[2]))
        R = sqrt(1 + ((p.Lc / p.L1) * tδ)^2)
        b2 = abs(p.Lc2) / p.L3

        bθ = v
        bϕ1 = (v / p.L2) * R
        b5ϕ1 = v * (R / p.L3 + (1 + b2) * R / p.L2)
        b5ϕ2 = v * (R / p.L3 + b2 * R / p.L2)

        return @SMatrix [
            0.0 0.0 bθ 0.0 0.0
            0.0 0.0 bθ 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 bϕ1 0.0
            0.0 0.0 0.0 b5ϕ1 b5ϕ2
        ]
    end
end

# Optional: scalar bound on ||A|| (crude but monotone).
function bound_norm_jacobian(p::Params = Params())
    return u -> begin
        v = abs(u[1])
        tδ = abs(tan(u[2]))
        R = sqrt(1 + ((p.Lc / p.L1) * tδ)^2)
        b2 = abs(p.Lc2) / p.L3
        return v +
               (v / p.L2) * R +
               v * (R / p.L3 + (1 + b2) * R / p.L2) +
               v * (R / p.L3 + b2 * R / p.L2)
    end
end

# Optional: bound on the Hessian tensor norm. The second ϕ-derivatives of f4
# and f5 are the same sin/cos combinations as the first (differentiation flips
# sin ↔ cos), so their harmonic amplitudes repeat — crude but monotone sum.
function bound_norm_hessian_tensor(p::Params = Params())
    return u -> begin
        v = abs(u[1])
        tδ = abs(tan(u[2]))
        R = sqrt(1 + ((p.Lc / p.L1) * tδ)^2)
        b2 = abs(p.Lc2) / p.L3
        return (v / p.L2) * R + v * (R / p.L3 + (1 + b2) * R / p.L2)
    end
end

# ----------------------------
# Dionysos system wrapper
# ----------------------------
function system(
    _X_;
    _U_ = LazySets.Hyperrectangle(; low = SVector(-3.0, -0.7), high = SVector(3.0, 0.7)),
    params::Params = Params(),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        LazySets.dim(_X_),
        LazySets.dim(_U_),
        _X_,
        _U_,
    )
end

function with_phi_limits(
    _X_::LazySets.AbstractHyperrectangle;
    phi1_max = 0.7,
    phi2_max = 0.7,
)
    lb = SVector(
        LazySets.low(_X_, 1),
        LazySets.low(_X_, 2),
        LazySets.low(_X_, 3),
        -phi1_max,
        -phi2_max,
    )
    ub = SVector(
        LazySets.high(_X_, 1),
        LazySets.high(_X_, 2),
        LazySets.high(_X_, 3),
        phi1_max,
        phi2_max,
    )
    return LazySets.Hyperrectangle(; low = lb, high = ub)
end

function extrude_xy_obstacle_to_5d(ob2d, _X_)
    lb = SVector(
        LazySets.low(ob2d, 1),
        LazySets.low(ob2d, 2),
        LazySets.low(_X_, 3),
        LazySets.low(_X_, 4),
        LazySets.low(_X_, 5),
    )
    ub = SVector(
        LazySets.high(ob2d, 1),
        LazySets.high(ob2d, 2),
        LazySets.high(_X_, 3),
        LazySets.high(_X_, 4),
        LazySets.high(_X_, 5),
    )
    return LazySets.Hyperrectangle(; low = lb, high = ub)
end

function with_xy_obstacles(
    _X_::LazySets.AbstractHyperrectangle;
    obstacles2d = corridor_obstacles(),
)
    obs5d = [extrude_xy_obstacle_to_5d(ob, _X_) for ob in obstacles2d]
    return UT.set_minus(_X_, UT.set_union(obs5d))
end

# ----------------------------
# Forward corridor-parking benchmark (replaces the PR's empty placeholder,
# which had no obstacles and was never run). Two walls leave a dead-end
# corridor (y ∈ [4.5, 7.5] for x ∈ [12, 30]) opening west, sized for the
# ~10 m rig (nose 1.25·L1 = 2.75 m ahead of the state point, trailer-2 axle
# Lc + L2 + Lc2 + L3 = 6.8 m behind). The rig starts in the open area heading
# east and drives nose-first into the corridor — the ϕ-STABLE direction, so
# both trailers follow instead of jackknifing.
# `margin` inflates the walls uniformly: the collision model is the tractor
# rear-axle POINT while the drawn bodies are ~1 m wide — a body-clearance
# margin keeps the certified point trajectory visually clear of the true walls.
_inflate(ob, margin) = LazySets.Hyperrectangle(;
    low = SVector(LazySets.low(ob, 1) - margin, LazySets.low(ob, 2) - margin),
    high = SVector(LazySets.high(ob, 1) + margin, LazySets.high(ob, 2) + margin),
)

corridor_obstacles(; margin = 0.0) = [
    _inflate(
        LazySets.Hyperrectangle(; low = SVector(12.0, 0.0), high = SVector(30.0, 4.5)),
        margin,
    ),
    _inflate(
        LazySets.Hyperrectangle(; low = SVector(12.0, 7.5), high = SVector(30.0, 12.0)),
        margin,
    ),
]

"""
    problem(; params, obstacles2d, phi_max) -> OptimalControlProblem

The two-trailer forward corridor-parking problem: domain `[0,30] × [0,12]`
minus the corridor walls, both hitch angles limited to `±phi_max`, inputs
`v ∈ [-3, 3]`, `δ ∈ [-0.7, 0.7]`. Initial states in the open area heading
east; target deep in the corridor (`[23,25] × [5.5,6.5]`, θ and both hitch
angles within ±5°) — parked poses whose 10 m body stays inside the corridor:
nose ≤ 27.75, trailer-2 axle ≥ 16.2.
"""
function problem(;
    params::Params = Params(),
    obstacles2d = corridor_obstacles(),
    phi_max = deg2rad(50.0),
    transition_cost = nothing,
    state_cost = nothing,
)
    box = LazySets.Hyperrectangle(;
        low = SVector(0.0, 0.0, -pi, -phi_max, -phi_max),
        high = SVector(30.0, 12.0, pi, phi_max, phi_max),
    )
    _X_ = with_xy_obstacles(box; obstacles2d = obstacles2d)

    _U_ = LazySets.Hyperrectangle(; low = SVector(-3.0, -0.7), high = SVector(3.0, 0.7))

    _I_ = LazySets.Hyperrectangle(;
        low = SVector(1.0, 1.0, -0.3, -0.2, -0.2),
        high = SVector(4.0, 3.0, 0.3, 0.2, 0.2),
    )
    _T_ = LazySets.Hyperrectangle(;
        low = SVector(23.0, 5.5, -deg2rad(5.0), -deg2rad(5.0), -deg2rad(5.0)),
        high = SVector(25.0, 6.5, deg2rad(5.0), deg2rad(5.0), deg2rad(5.0)),
    )

    sys = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        5,
        2,
        _X_,
        _U_,
    )
    return PR.OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost)
end

################################################
########## Visualization tools #################
################################################

Base.@kwdef struct DrawParams{T}
    tractor_length::T
    tractor_width::T
    trailer1_length::T
    trailer2_length::T
    trailer_width::T
    axle_halfwidth::T
end
function DrawParams(
    p::Params{T};
    tractor_length = 1.25 * p.L1,
    tractor_width = 0.45 * p.L1,
    trailer1_length = 1.00 * p.L2,
    trailer2_length = 1.00 * p.L3,
    trailer_width = 0.40 * p.L1,
    axle_halfwidth = 0.25 * p.L1,
) where {T}
    return DrawParams{T}(
        tractor_length,
        tractor_width,
        trailer1_length,
        trailer2_length,
        trailer_width,
        axle_halfwidth,
    )
end

rot2(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

function rect_poly(c::SVector{2, Float64}, θ, L, W)
    R = rot2(θ)
    hl, hw = L / 2, W / 2
    pts = (
        c + R * SVector(hl, hw),
        c + R * SVector(hl, -hw),
        c + R * SVector(-hl, -hw),
        c + R * SVector(-hl, hw),
        c + R * SVector(hl, hw),
    )
    return [p[1] for p in pts], [p[2] for p in pts]
end

# Forward-kinematics chain of the three bodies.
function vehicle_keypoints(p::Params, x)
    θ1 = x[3]
    θ2 = θ1 + x[4]
    θ3 = θ2 + x[5]
    P_rear = SVector(x[1], x[2])
    P_front = P_rear + p.L1 * SVector(cos(θ1), sin(θ1))
    H1 = P_rear - p.Lc * SVector(cos(θ1), sin(θ1))
    P_tr1 = H1 - p.L2 * SVector(cos(θ2), sin(θ2))
    H2 = P_tr1 - p.Lc2 * SVector(cos(θ2), sin(θ2))
    P_tr2 = H2 - p.L3 * SVector(cos(θ3), sin(θ3))
    return P_rear, P_front, H1, P_tr1, H2, P_tr2, θ1, θ2, θ3
end

function draw_articulated!(plt, p::Params, dp::DrawParams, x, u)
    P_rear, P_front, H1, P_tr1, H2, P_tr2, θ1, θ2, θ3 = vehicle_keypoints(p, x)

    tractor_center = P_rear + 0.5 * dp.tractor_length * SVector(cos(θ1), sin(θ1))
    tx, ty = rect_poly(tractor_center, θ1, dp.tractor_length, dp.tractor_width)
    plot!(plt, tx, ty; lw = 1, fill = (true, 0.10), label = false)

    tr1_center = P_tr1 + 0.5 * dp.trailer1_length * SVector(cos(θ2), sin(θ2))
    r1x, r1y = rect_poly(tr1_center, θ2, dp.trailer1_length, dp.trailer_width)
    plot!(plt, r1x, r1y; lw = 1, fill = (true, 0.08), label = false)

    tr2_center = P_tr2 + 0.5 * dp.trailer2_length * SVector(cos(θ3), sin(θ3))
    r2x, r2y = rect_poly(tr2_center, θ3, dp.trailer2_length, dp.trailer_width)
    plot!(plt, r2x, r2y; lw = 1, fill = (true, 0.08), label = false)

    # Hitch links.
    plot!(plt, [P_rear[1], H1[1]], [P_rear[2], H1[2]]; lw = 2, label = false)
    plot!(plt, [P_tr1[1], H2[1]], [P_tr1[2], H2[2]]; lw = 2, label = false)

    # Axles.
    for (P, θ) in ((P_rear, θ1), (P_front, θ1), (P_tr1, θ2), (P_tr2, θ3))
        n = SVector(-sin(θ), cos(θ))
        a = P + dp.axle_halfwidth * n
        b = P - dp.axle_halfwidth * n
        plot!(plt, [a[1], b[1]], [a[2], b[2]]; lw = 2, label = false)
    end

    return plt
end

function plot_xy_obstacles!(plt, obs2d; alpha = 0.8, color = :black, kwargs...)
    for ob in obs2d
        x1l, y1l = LazySets.low(ob, 1), LazySets.low(ob, 2)
        x1u, y1u = LazySets.high(ob, 1), LazySets.high(ob, 2)
        xs = [x1l, x1u, x1u, x1l, x1l]
        ys = [y1l, y1l, y1u, y1u, y1l]
        plot!(
            plt,
            xs,
            ys;
            lw = 1,
            color = color,
            fill = (true, alpha),
            label = false,
            kwargs...,
        )
    end
    return plt
end

function system_plot!(;
    params::Params = Params(),
    draw_params::DrawParams = DrawParams(params),
    xlims = nothing,
    ylims = nothing,
    obstacles2d = nothing,
)
    return function (fig, x, u)
        obstacles2d !== nothing && plot_xy_obstacles!(fig, obstacles2d)
        draw_articulated!(fig, params, draw_params, x, u)
        xlims !== nothing && xlims!(fig, xlims...)
        ylims !== nothing && ylims!(fig, ylims...)
        return fig
    end
end

end # module
