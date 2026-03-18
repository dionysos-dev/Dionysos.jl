module ArticulatedVehicle2Trailers

using StaticArrays
using MathematicalSystems
using Dionysos
using Plots
import Symbolics
import IntervalArithmetic as IA

const UT = Dionysos.Utils
const ST = Dionysos.System
const PB = Dionysos.Problem
const SY = Dionysos.Symbolic

# ----------------------------
# Parameters 
# ----------------------------
Base.@kwdef struct Params{T}
    L1::T = 2.0      # tractor wheelbase
    L2::T = 6.0      # hitch 1 -> trailer 1 axle distance
    L3::T = 6.0      # hitch 2 -> trailer 2 axle distance
    Lc::T = 1.0      # off-axle hitch length behind the tractor
    Lc2::T = 0.0     # off-axle hitch length behind trailer 1
end

# ----------------------------
# Dynamics: ẋ = f(x,u)
# x = [x1, x2, θ1, ϕ1, ϕ2]
# u = [v1, σ̃] with σ̃ = tan(δ)
# ----------------------------
function dynamic(p::Params = Params())
    return (x, u) -> begin
        v = u[1]
        σ̃ = u[2]
        θ1 = x[3]
        ϕ1 = x[4]
        ϕ2 = x[5]

        return SVector{5}(
            v * cos(θ1),
            v * sin(θ1),
            (v / p.L1) * σ̃,
            -(v / (p.L1 * p.L2)) * (p.L1 * sin(ϕ1) + (p.Lc * cos(ϕ1) + p.L2) * σ̃), 
            -(v/p.L3) * (cos(ϕ1) - (p.Lc/p.L1) * σ̃ * sin(ϕ1)) * sin(ϕ2) + (1 + (p.Lc2/p.L3) * cos(ϕ2)) * ( (v/p.L2) * sin(ϕ1) + ((v * p.Lc)/(p.L1 * p.L2)) * σ̃ * cos(ϕ1))
        )
    end
end

# ----------------------------
# Jacobian A(x,u) = ∂f/∂x
# ----------------------------
function jacobian(p::Params = Params())
    return (x, u) -> begin
        v = u[1]
        σ̃ = u[2]
        θ1 = x[3]
        ϕ1 = x[4]
        ϕ2 = x[5]

        A = cos(ϕ1) - (p.Lc / p.L1) * σ̃ * sin(ϕ1)
        B = (1 / p.L2) * sin(ϕ1) + (p.Lc / (p.L1 * p.L2)) * σ̃ * cos(ϕ1)

        Aϕ = -sin(ϕ1) - (p.Lc / p.L1) * σ̃ * cos(ϕ1)
        Bϕ = (1 / p.L2) * cos(ϕ1) - (p.Lc / (p.L1 * p.L2)) * σ̃ * sin(ϕ1)

        d4dϕ1 = -(v / (p.L1 * p.L2)) * (p.L1 * cos(ϕ1) - p.Lc * σ̃ * sin(ϕ1))
        d5dϕ1 = -(v / p.L3) * Aϕ * sin(ϕ2) + v * (1 + (p.Lc2 / p.L3) * cos(ϕ2)) * Bϕ
        d5dϕ2 = -(v / p.L3) * A * cos(ϕ2) - v * (p.Lc2 / p.L3) * B * sin(ϕ2)

        return @SMatrix [
            0.0 0.0 -v * sin(θ1) 0.0 0.0
            0.0 0.0  v * cos(θ1) 0.0 0.0
            0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 d4dϕ1 0.0
            0.0 0.0 0.0 d5dϕ1 d5dϕ2
        ]
    end
end

# ----------------------------
# A cheap bound matrix Jb(u) such that |A(x,u)| ≤ Jb(u) elementwise
# (useful for abstraction growth bounds)
# ----------------------------
# function jacobian_bound(p::Params = Params())
#     return u -> begin
#         v = abs(u[1])
#         σ̃ = abs(u[2])

#         bθ = v
#         bϕ1 = v / (p.L1 * p.L2) * (p.L1 + abs(p.Lc) * σ̃)
#         bϕ2_ϕ1 = v * (
#             (1 + abs(p.Lc) * σ̃ / p.L1) / p.L3 +
#             (1 + abs(p.Lc2) / p.L3) * (1 / p.L2 + abs(p.Lc) * σ̃ / (p.L1 * p.L2))
#         )
#         bϕ2_ϕ2 = v / p.L3 * (
#             1 + abs(p.Lc) * σ̃ / p.L1 +
#             abs(p.Lc2) * (1 / p.L2 + abs(p.Lc) * σ̃ / (p.L1 * p.L2))
#         )

#         return @SMatrix [
#             0.0 0.0 bθ 0.0 0.0
#             0.0 0.0 bθ 0.0 0.0
#             0.0 0.0 0.0 0.0 0.0
#             0.0 0.0 0.0 bϕ1 0.0
#             0.0 0.0 0.0 bϕ2_ϕ1 bϕ2_ϕ2
#         ]
#     end
# end

function jacobian_bound(p::Params = Params())
    return u -> begin
        v = abs(u[1])
        σ = u[2]

        a = (p.Lc / p.L1) * σ
        R = sqrt(1 + a^2)
        b = p.Lc2 / p.L3

        bθ = v
        bϕ1 = (v / p.L2) * R

        bϕ2_ϕ1 = v * (
            R / p.L3 +
            (1 + abs(b)) * R / p.L2
        )

        bϕ2_ϕ2 = v * (
            R / p.L3 +
            abs(b) * R / p.L2
        )

        return @SMatrix [
            0.0  0.0  bθ       0.0        0.0
            0.0  0.0  bθ       0.0        0.0
            0.0  0.0  0.0      0.0        0.0
            0.0  0.0  0.0      bϕ1        0.0
            0.0  0.0  0.0      bϕ2_ϕ1     bϕ2_ϕ2
        ]
    end
end

function bound_norm_jacobian(p::Params = Params())
    return u -> begin
        v = abs(u[1])
        σ̃ = abs(u[2])
        bϕ1 = v / (p.L1 * p.L2) * (p.L1 + abs(p.Lc) * σ̃)
        bϕ2_ϕ1 = v * (
            (1 + abs(p.Lc) * σ̃ / p.L1) / p.L3 +
            (1 + abs(p.Lc2) / p.L3) * (1 / p.L2 + abs(p.Lc) * σ̃ / (p.L1 * p.L2))
        )
        bϕ2_ϕ2 = v / p.L3 * (
            1 + abs(p.Lc) * σ̃ / p.L1 +
            abs(p.Lc2) * (1 / p.L2 + abs(p.Lc) * σ̃ / (p.L1 * p.L2))
        )
        # crude but monotone (you can replace by something tighter)
        return v + bϕ1 + bϕ2_ϕ1 + bϕ2_ϕ2
    end
end

function bound_norm_hessian_tensor(p::Params = Params())
    return u -> begin
        v = abs(u[1])
        σ̃ = abs(u[2])
        bϕ1 = v / (p.L1 * p.L2) * (p.L1 + abs(p.Lc) * σ̃)
        bϕ2_ϕ1 = v * (
            (1 + abs(p.Lc) * σ̃ / p.L1) / p.L3 +
            (1 + abs(p.Lc2) / p.L3) * (1 / p.L2 + abs(p.Lc) * σ̃ / (p.L1 * p.L2))
        )
        bϕ2_ϕ2 = v / p.L3 * (
            1 + abs(p.Lc) * σ̃ / p.L1 +
            abs(p.Lc2) * (1 / p.L2 + abs(p.Lc) * σ̃ / (p.L1 * p.L2))
        )
        return bϕ1 + bϕ2_ϕ1 + bϕ2_ϕ2
    end
end

# ----------------------------
# Dionysos system wrapper
# ----------------------------
function system(
    _X_;
    _U_ = UT.HyperRectangle(SVector(-1.0, -0.6), SVector(1.0, 0.6)),
    params::Params = Params(),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

"""
    symbolic_system(_X_; _U_, params, Ts, ΔX, ΔU, ΔW)

Version symbolique discrete pour la synthese LMI.
Le systeme black-box reste utilise pour la simulation nominale.

"""
function symbolic_system(
    _X_;
    _U_ = UT.HyperRectangle(SVector(-1.0, -0.6), SVector(1.0, 0.6)),
    params::Params = Params(),
    Ts::Float64 = 0.001,
    ΔX = IA.IntervalBox(IA.interval(-0.2, 0.2), 5),
    ΔU = IA.IntervalBox(IA.interval(-0.15, 0.15), 2),
    ΔW = IA.IntervalBox(IA.interval(0.0, 0.0), 1),
    rk4_num_substeps::Int = 1,
    obstacles = Any[],
)
    Symbolics.@variables x1 x2 θ1 ϕ1 ϕ2 v1 σ̃ w1 T
    x = [x1; x2; θ1; ϕ1; ϕ2]
    u = [v1; σ̃]
    w = [w1]

    f_cont_fun = dynamic(params)
    f_cont_expr(xloc, uloc) = collect(f_cont_fun(xloc, uloc))


    # Discretisation symbolique choisie pour la linearisation LMI.
    f_disc = ST.runge_kutta4(f_cont_expr, x, u, T, rk4_num_substeps)

    fsymbolicT = eval(ST.build_function(f_disc, x, u, w, T)[1])
    fsymbolic = Symbolics.substitute(f_disc, Dict(T => Ts))

    # No additive noise in this model, but LMI API expects a noise format.
    Wset = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    Uformat = UT.format_input_set(_U_)
    Wformat = UT.format_noise_set(Wset)

    f_cont_fun = dynamic(params)
    function f_eval(xv, uv, _wv)
        xsv = SVector{5, Float64}(xv)
        usv = SVector{2, Float64}(uv)
        xnext = ST.runge_kutta4(f_cont_fun, xsv, usv, Ts, rk4_num_substeps)
        return collect(xnext)
    end
    function f_backward_eval(xv, uv, _wv)
        xsv = SVector{5, Float64}(xv)
        usv = SVector{2, Float64}(uv)
        xprev = ST.runge_kutta4(f_cont_fun, xsv, usv, -Ts, rk4_num_substeps)
        return collect(xprev)
    end

    return ST.SymbolicSystem(
        fsymbolicT,
        fsymbolic,
        Ts,
        length(x),
        length(u),
        length(w),
        x,
        u,
        w,
        ΔX,
        ΔU,
        ΔW,
        _X_,
        _U_,
        Wset,
        obstacles,
        f_eval,
        f_backward_eval,
        Uformat,
        Wformat,
    )
end

function with_phi_limits(_X_::UT.HyperRectangle; phi1_max = 0.7, phi2_max = 0.7)
    lb = SVector(_X_.lb[1], _X_.lb[2], _X_.lb[3], -phi1_max, -phi2_max)
    ub = SVector(_X_.ub[1], _X_.ub[2], _X_.ub[3], phi1_max, phi2_max)
    return UT.HyperRectangle(lb, ub)
end

function with_phi_limit(_X_::UT.HyperRectangle; phi_max = 0.7)
    return with_phi_limits(_X_; phi1_max = phi_max, phi2_max = phi_max)
end

function extrude_xy_obstacle_to_5d(ob2d, _X_)
    # ob2d has lb=[x,y], ub=[x,y]
    lb = SVector(ob2d.lb[1], ob2d.lb[2], _X_.lb[3], _X_.lb[4], _X_.lb[5])
    ub = SVector(ob2d.ub[1], ob2d.ub[2], _X_.ub[3], _X_.ub[4], _X_.ub[5])
    return UT.HyperRectangle(lb, ub)
end

function with_xy_obstacles(_X_::UT.HyperRectangle; obstacles2d = Any[])
    isempty(obstacles2d) && return _X_
    obs5d = [extrude_xy_obstacle_to_5d(ob, _X_) for ob in obstacles2d]
    return UT.LazySetMinus(_X_, UT.LazySetUnion(obs5d))
end

# ----------------------------
# Benchmark "problem" factory
# (you can tune X/I/T to your benchmark)
# ----------------------------
function problem(;
    params::Params = Params(),
    transition_cost = nothing,
    state_cost = nothing,
    terminal_cost = PB.Infinity(),
)

    # Example domains (edit):
    _X_ = UT.HyperRectangle(
        SVector(-20.0, -20.0, -pi, -pi / 2, -pi / 2),
        SVector(20.0, 20.0, pi, pi / 2, pi / 2),
    )

    _I_ = UT.HyperRectangle(
        SVector(-10.0, -10.0, -0.2, 0.0, 0.0),
        SVector(-9.0, -9.0, 0.2, 0.0, 0.0),
    )

    _T_ = UT.HyperRectangle(
        SVector(9.0, 9.0, -0.2, -0.2, -0.2),
        SVector(10.0, 10.0, 0.2, 0.2, 0.2),
    )

    sys = system(_X_; params = params)

    # If you want pure reachability: state_cost = nothing, transition_cost = nothing.
    return PB.OptimalControlProblem(
        sys,
        _I_,
        _T_,
        state_cost,
        transition_cost,
        terminal_cost,
    )
end

################################################
############ Simple Controllers ################
################################################

function get_constant_controller(u_const)
    return ST.ConstantController(u_const)
end

function get_goal_seeking_controller(xg, yg; v = 1.0, δmax = 0.5, k = 1.2)
    f = x -> begin
        x1, x2, θ, ϕ1, ϕ2 = x
        desired = atan(yg - x2, xg - x1)
        e = mod(desired - θ + pi, 2pi) - pi
        δ = clamp(k * e, -δmax, δmax)
        return SVector(v, tan(δ))
    end
    return ST.BlackBoxContinuousController(f)
end

################################################
########## Visualization tools #################
################################################
Base.@kwdef struct DrawParams{T}
    tractor_length::T
    tractor_width::T
    trailer1_length::T
    trailer1_width::T
    trailer2_length::T
    trailer2_width::T
    wheel_length::T
    wheel_width::T
    axle_halfwidth::T
end

function DrawParams(
    p::Params{T};
    tractor_length = 1.25 * p.L1,      # body a bit longer than wheelbase
    tractor_width = 0.45 * p.L1,
    trailer1_length = 1.00 * p.L2,     # trailer body ≈ trailer length
    trailer1_width = 0.40 * p.L1,
    trailer2_length = 1.00 * p.L3,
    trailer2_width = 0.40 * p.L1,
    wheel_length = 0.22 * p.L1,
    wheel_width = 0.08 * p.L1,
    axle_halfwidth = 0.25 * p.L1,
) where {T}
    return DrawParams{T}(
        tractor_length,
        tractor_width,
        trailer1_length,
        trailer1_width,
        trailer2_length,
        trailer2_width,
        wheel_length,
        wheel_width,
        axle_halfwidth,
    )
end

rot2(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

# rectangle centered at c, oriented by θ, size (L, W)
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
    xs = [p[1] for p in pts]
    ys = [p[2] for p in pts]
    return xs, ys
end

function vehicle_keypoints(p, x, dp::DrawParams)
    x1, x2, θ1, ϕ1, ϕ2 = x
    θ2 = θ1 + ϕ1
    θ3 = θ2 + ϕ2

    P1 = SVector(x1, x2)                                   # tractor rear axle (state position)
    Pfront = P1 + p.L1 * SVector(cos(θ1), sin(θ1))         # tractor front axle

    H1 = P1 - p.Lc * SVector(cos(θ1), sin(θ1))             # hitch behind rear axle
    P2 = H1 - p.L2 * SVector(cos(θ2), sin(θ2))             # trailer 1 axle center

    H2 = P2 - p.Lc2 * SVector(cos(θ2), sin(θ2))            # hitch behind trailer 1 axle
    P3 = H2 - p.L3 * SVector(cos(θ3), sin(θ3))             # trailer 2 axle center

    return P1, Pfront, H1, P2, H2, P3, θ1, θ2, θ3
end

function draw_articulated!(
    plt,
    p,
    dp::DrawParams,
    x,
    u;
    show_axes = true,
    show_heading = true,
    show_phi_arc = true,
)
    P1, Pfront, H1, P2, H2, P3, θ1, θ2, θ3 = vehicle_keypoints(p, x, dp)
    σ̃ = u[2]
    δ = atan(σ̃)

    # --- Bodies (rectangles) ---
    tractor_center = P1 + 0.5 * dp.tractor_length * SVector(cos(θ1), sin(θ1))
    tx, ty = rect_poly(tractor_center, θ1, dp.tractor_length, dp.tractor_width)
    plot!(plt, tx, ty; lw = 1, fill = (true, 0.10), label = false)

    trailer1_center = P2 + 0.5 * dp.trailer1_length * SVector(cos(θ2), sin(θ2))
    rx1, ry1 = rect_poly(trailer1_center, θ2, dp.trailer1_length, dp.trailer1_width)
    plot!(plt, rx1, ry1; lw = 1, fill = (true, 0.08), label = false)

    trailer2_center = P3 + 0.5 * dp.trailer2_length * SVector(cos(θ3), sin(θ3))
    rx2, ry2 = rect_poly(trailer2_center, θ3, dp.trailer2_length, dp.trailer2_width)
    plot!(plt, rx2, ry2; lw = 1, fill = (true, 0.06), label = false)

    # --- Axles (lines) ---
    if show_axes
        n1 = SVector(-sin(θ1), cos(θ1))
        A1a = P1 + dp.axle_halfwidth * n1
        A1b = P1 - dp.axle_halfwidth * n1
        plot!(plt, [A1a[1], A1b[1]], [A1a[2], A1b[2]]; lw = 2, label = false)

        A2a = Pfront + dp.axle_halfwidth * n1
        A2b = Pfront - dp.axle_halfwidth * n1
        plot!(plt, [A2a[1], A2b[1]], [A2a[2], A2b[2]]; lw = 2, label = false)

        n2 = SVector(-sin(θ2), cos(θ2))
        T1a = P2 + dp.axle_halfwidth * n2
        T1b = P2 - dp.axle_halfwidth * n2
        plot!(plt, [T1a[1], T1b[1]], [T1a[2], T1b[2]]; lw = 2, label = false)

        n3 = SVector(-sin(θ3), cos(θ3))
        T2a = P3 + dp.axle_halfwidth * n3
        T2b = P3 - dp.axle_halfwidth * n3
        plot!(plt, [T2a[1], T2b[1]], [T2a[2], T2b[2]]; lw = 2, label = false)
    end

    # --- Wheels (rectangles) ---
    n1 = SVector(-sin(θ1), cos(θ1))
    for s in (-1.0, 1.0)
        c = P1 + s * dp.axle_halfwidth * n1
        wx, wy = rect_poly(c, θ1, dp.wheel_length, dp.wheel_width)
        plot!(plt, wx, wy; lw = 1, fill = (true, 0.25), label = false)
    end

    θw = θ1 + δ
    for s in (-1.0, 1.0)
        c = Pfront + s * dp.axle_halfwidth * n1
        wx, wy = rect_poly(c, θw, dp.wheel_length, dp.wheel_width)
        plot!(plt, wx, wy; lw = 1, fill = (true, 0.25), label = false)
    end

    n2 = SVector(-sin(θ2), cos(θ2))
    for s in (-1.0, 1.0)
        c = P2 + s * dp.axle_halfwidth * n2
        wx, wy = rect_poly(c, θ2, dp.wheel_length, dp.wheel_width)
        plot!(plt, wx, wy; lw = 1, fill = (true, 0.25), label = false)
    end

    n3 = SVector(-sin(θ3), cos(θ3))
    for s in (-1.0, 1.0)
        c = P3 + s * dp.axle_halfwidth * n3
        wx, wy = rect_poly(c, θ3, dp.wheel_length, dp.wheel_width)
        plot!(plt, wx, wy; lw = 1, fill = (true, 0.25), label = false)
    end

    # --- Hitch links ---
    plot!(plt, [P1[1], H1[1]], [P1[2], H1[2]]; lw = 2, label = false)
    plot!(plt, [P2[1], H2[1]], [P2[2], H2[2]]; lw = 2, label = false)

    # --- Heading arrows ---
    if show_heading
        a1 = P1 + 1.2 * SVector(cos(θ1), sin(θ1))
        plot!(plt, [P1[1], a1[1]], [P1[2], a1[2]]; lw = 2, label = false)

        a2 = P2 + 1.2 * SVector(cos(θ2), sin(θ2))
        plot!(plt, [P2[1], a2[1]], [P2[2], a2[2]]; lw = 2, label = false)

        a3 = P3 + 1.2 * SVector(cos(θ3), sin(θ3))
        plot!(plt, [P3[1], a3[1]], [P3[2], a3[2]]; lw = 2, label = false)

        sdir = Pfront + 1.0 * SVector(cos(θw), sin(θw))
        plot!(
            plt,
            [Pfront[1], sdir[1]],
            [Pfront[2], sdir[2]];
            lw = 2,
            ls = :dash,
            label = false,
        )
    end

    # --- φ arcs (simple polyline arcs around hitches) ---
    if show_phi_arc
        ϕ1 = x[4]
        ϕ2 = x[5]

        r1 = 1.0
        ts1 = range(θ1, θ1 + ϕ1; length = 25)
        ax1 = [H1[1] + r1 * cos(t) for t in ts1]
        ay1 = [H1[2] + r1 * sin(t) for t in ts1]
        plot!(plt, ax1, ay1; lw = 2, ls = :dot, label = false)

        r2 = 0.8
        ts2 = range(θ2, θ2 + ϕ2; length = 25)
        ax2 = [H2[1] + r2 * cos(t) for t in ts2]
        ay2 = [H2[2] + r2 * sin(t) for t in ts2]
        plot!(plt, ax2, ay2; lw = 2, ls = :dot, label = false)
    end

    return plt
end

function plot_xy_obstacles!(plt, obs2d; alpha = 0.25, kwargs...)
    for ob in obs2d
        x1l, y1l = ob.lb[1], ob.lb[2]
        x1u, y1u = ob.ub[1], ob.ub[2]
        xs = [x1l, x1u, x1u, x1l, x1l]
        ys = [y1l, y1l, y1u, y1u, y1l]
        plot!(plt, xs, ys; lw = 1, fill = (true, alpha), label = false, kwargs...)
    end
    return plt
end

function live_vehicle_progression(
    p,
    dp,
    x_traj::ST.Trajectory,
    u_traj::ST.Trajectory,
    xl,
    yl;
    domain = nothing,
    obstacles2d = nothing,
    every = 1,
    dt = 0.05,
    giffile::Union{Nothing, String} = nothing,
    fps::Int = 20,
    title::Union{Nothing, String} = nothing,
)
    states = x_traj.seq
    inputs = u_traj.seq

    xs = [x[1] for x in states]
    ys = [x[2] for x in states]

    plot_title = title === nothing ? "" : title

    # --- GIF MODE ---
    if giffile !== nothing
        anim = @animate for k in 1:every:length(states)
            plt = plot(;
                aspect_ratio = :equal,
                xlims = xl,
                ylims = yl,
                legend = false,
                size = (700, 700),
                title = plot_title,
            )
            if domain !== nothing
                plot!(plt, domain; color = :grey, opacity = 0.1)
            end
            if obstacles2d !== nothing
                plot_xy_obstacles!(plt, obstacles2d; color = :black)
            end
            plot!(plt, xs, ys; lw = 1)
            uk = (k <= length(inputs)) ? inputs[k] : inputs[end]
            draw_articulated!(plt, p, dp, states[k], uk)
        end

        gif(anim, giffile; fps = fps)
        return anim
    end

    # --- LIVE MODE ---
    for k in 1:every:length(states)
        plt = plot(;
            aspect_ratio = :equal,
            xlims = xl,
            ylims = yl,
            legend = false,
            size = (700, 700),
            title = plot_title,
        )
        if domain !== nothing
            plot!(plt, domain; color = :grey, opacity = 0.1)
        end
        if obstacles2d !== nothing
            plot_xy_obstacles!(plt, obstacles2d; color = :black)
        end
        plot!(plt, xs, ys; lw = 1)
        uk = (k <= length(inputs)) ? inputs[k] : inputs[end]
        draw_articulated!(plt, p, dp, states[k], uk)
        display(plt)
        sleep(dt)
    end

    return nothing
end

end # module
