module ArticulatedVehicle

using StaticArrays
using MathematicalSystems
using Dionysos
using Plots

const UT = Dionysos.Utils
const ST = Dionysos.System
const PB = Dionysos.Problem

# ----------------------------
# Parameters (edit as you want)
# ----------------------------
Base.@kwdef struct Params{T}
    L1::T = 2.0      # tractor wheelbase
    L2::T = 6.0      # trailer length
    Lc::T = 1.0      # off-axle hitch length (called "c" in the paper)
end

# ----------------------------
# Dynamics: ẋ = f(x,u)
# x = [x1, x2, θ1, ϕ]
# u = [v1, δ]
# ----------------------------
function dynamic(p::Params = Params())
    return (x, u) -> begin
        v  = u[1]
        δ  = u[2]
        θ  = x[3]
        ϕ  = x[4]

        return SVector{4}(
            v * cos(θ),
            v * sin(θ),
            (v / p.L1) * tan(δ),
            -(v / (p.L1 * p.L2)) * (p.L1 * sin(ϕ) + (p.Lc * cos(ϕ) + p.L2) * tan(δ)),
        )
    end
end

# ----------------------------
# Jacobian A(x,u) = ∂f/∂x
# ----------------------------
function jacobian(p::Params = Params())
    return (x, u) -> begin
        v  = u[1]
        θ  = x[3]
        ϕ  = x[4]
        δ  = u[2]
        tδ = tan(δ)

        d4dϕ = -(v / (p.L1 * p.L2)) * (p.L1 * cos(ϕ) - p.Lc * sin(ϕ) * tδ)
        return SMatrix{4,4}(
            0.0, 0.0, -v*sin(θ), 0.0,
            0.0, 0.0,  v*cos(θ), 0.0,
            0.0, 0.0,      0.0,  0.0,
            0.0, 0.0,      0.0,  d4dϕ,
        )
    end
end

# ----------------------------
# A cheap bound matrix Jb(u) such that |A(x,u)| ≤ Jb(u) elementwise
# (useful for abstraction growth bounds)
# ----------------------------
function jacobian_bound(p::Params = Params())
    return u -> begin
        v  = abs(u[1])
        δ  = u[2]
        tδ = abs(tan(δ))

        # Bounds:
        # |∂f1/∂θ| ≤ |v|, |∂f2/∂θ| ≤ |v|
        # |∂f4/∂ϕ| ≤ |v|/(L1 L2)*(L1*|cosϕ| + |Lc|*|sinϕ|*|tanδ|)
        #          ≤ |v|/(L1 L2)*(L1 + |Lc|*|tanδ|)
        bθ = v
        bϕ = v/(p.L1*p.L2) * (p.L1 + abs(p.Lc)*tδ)

        return SMatrix{4,4}(
            0.0, 0.0, bθ,  0.0,
            0.0, 0.0, bθ,  0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, bϕ,
        )
    end
end

# Optional: scalar bound on ||A|| (simple conservative one)
function bound_norm_jacobian(p::Params = Params())
    return u -> begin
        v  = abs(u[1])
        tδ = abs(tan(u[2]))
        bϕ = v/(p.L1*p.L2) * (p.L1 + abs(p.Lc)*tδ)
        # crude but monotone (you can replace by something tighter)
        return v + bϕ
    end
end

# Optional: bound on Hessian tensor norm (only f4 has nonzero ϕϕ term)
function bound_norm_hessian_tensor(p::Params = Params())
    return u -> begin
        v  = abs(u[1])
        tδ = abs(tan(u[2]))
        # |∂²f4/∂ϕ²| ≤ |v|/(L1 L2)*(L1 + |Lc|*|tanδ|)
        return v/(p.L1*p.L2) * (p.L1 + abs(p.Lc)*tδ)
    end
end

# ----------------------------
# Dionysos system wrapper
# ----------------------------
function system(_X_;
    _U_ = UT.HyperRectangle(SVector(-1.0, -0.6), SVector(1.0, 0.6)),
    params::Params = Params()
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function with_phi_limit(_X_::UT.HyperRectangle; phi_max=0.7)
    lb = SVector(_X_.lb[1], _X_.lb[2], _X_.lb[3], -phi_max)
    ub = SVector(_X_.ub[1], _X_.ub[2], _X_.ub[3],  phi_max)
    return UT.HyperRectangle(lb, ub)
end

function extrude_xy_obstacle_to_4d(ob2d, _X_)
    # ob2d has lb=[x,y], ub=[x,y]
    lb = SVector(ob2d.lb[1], ob2d.lb[2], _X_.lb[3], _X_.lb[4])
    ub = SVector(ob2d.ub[1], ob2d.ub[2], _X_.ub[3], _X_.ub[4])
    return UT.HyperRectangle(lb, ub)
end
function with_xy_obstacles(_X_::UT.HyperRectangle; obstacles2d=xy_obstacles())
    obs4d = [extrude_xy_obstacle_to_4d(ob, _X_) for ob in obstacles2d]
    return UT.LazySetMinus(_X_, UT.LazyUnionSetArray(obs4d))
end

# ----------------------------
# Benchmark "problem" factory
# (you can tune X/I/T to your benchmark)
# ----------------------------
function problem(; params::Params = Params(),
                   transition_cost = nothing,
                   state_cost = nothing,
                   terminal_cost = PB.Infinity())

    # Example domains (edit):
    _X_ = UT.HyperRectangle(
        SVector(-20.0, -20.0, -pi, -pi/2),
        SVector( 20.0,  20.0,  pi,  pi/2),
    )

    _I_ = UT.HyperRectangle(
        SVector(-10.0, -10.0, -0.2, 0.0),
        SVector( -9.0,  -9.0,  0.2, 0.0),
    )

    _T_ = UT.HyperRectangle(
        SVector(  9.0,   9.0, -0.2, -0.2),
        SVector( 10.0,  10.0,  0.2,  0.2),
    )

    sys = system(_X_; params=params)

    # If you want pure reachability: state_cost = nothing, transition_cost = nothing.
    return PB.OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost, terminal_cost)
end

################################################
############ Simple Controllers ################
################################################

function get_constant_controller(u_const)
    return ST.ConstantController(u_const)
end

function get_goal_seeking_controller(xg, yg; v=1.0, δmax=0.5, k=1.2)
    f = x -> begin
        x1,x2,θ,ϕ = x
        desired = atan(yg - x2, xg - x1)
        e = mod(desired - θ + pi, 2pi) - pi
        δ = clamp(k*e, -δmax, δmax)
        return SVector(v, δ)
    end
    return ST.BlackBoxContinuousController(f)
end

################################################
########## Visualization tools #################
################################################
Base.@kwdef struct DrawParams{T}
    tractor_length::T
    tractor_width::T
    trailer_length::T
    trailer_width::T
    wheel_length::T
    wheel_width::T
    axle_halfwidth::T
end
function DrawParams(p::Params{T};
    tractor_length = 1.25*p.L1,          # body a bit longer than wheelbase
    tractor_width  = 0.45*p.L1,
    trailer_length = 1.00*p.L2,          # trailer body ≈ trailer length
    trailer_width  = 0.40*p.L1,          # often similar width as tractor
    wheel_length   = 0.22*p.L1,
    wheel_width    = 0.08*p.L1,
    axle_halfwidth = 0.25*p.L1,
) where {T}
    return DrawParams{T}(
        tractor_length,
        tractor_width,
        trailer_length,
        trailer_width,
        wheel_length,
        wheel_width,
        axle_halfwidth,
    )
end

rot2(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]

# rectangle centered at c, oriented by θ, size (L, W)
function rect_poly(c::SVector{2,Float64}, θ, L, W)
    R = rot2(θ)
    hl, hw = L/2, W/2
    pts = (
        c + R*SVector( hl,  hw),
        c + R*SVector( hl, -hw),
        c + R*SVector(-hl, -hw),
        c + R*SVector(-hl,  hw),
        c + R*SVector( hl,  hw),
    )
    xs = [p[1] for p in pts]
    ys = [p[2] for p in pts]
    return xs, ys
end

function vehicle_keypoints(p, x, dp::DrawParams)
    x1, x2, θ1, ϕ = x
    θ2 = θ1 + ϕ

    P_rear = SVector(x1, x2)                          # tractor rear axle (state position)
    P_front = P_rear + p.L1 * SVector(cos(θ1), sin(θ1))  # tractor front axle (wheelbase L1)

    H = P_rear - p.Lc * SVector(cos(θ1), sin(θ1))      # hitch behind rear axle
    P_tr = H - p.L2 * SVector(cos(θ2), sin(θ2))        # trailer axle center

    return P_rear, P_front, H, P_tr, θ1, θ2
end

function draw_articulated!(
    plt,
    p, dp::DrawParams,
    x, u;
    show_axes=true,
    show_heading=true,
    show_phi_arc=true,
)
    P_rear, P_front, H, P_tr, θ1, θ2 = vehicle_keypoints(p, x, dp)
    δ = u[2]

    # --- Bodies (rectangles) ---
    # Tractor body: center somewhere between rear and front (rough)
    tractor_center = P_rear + 0.5*dp.tractor_length * SVector(cos(θ1), sin(θ1))
    tx, ty = rect_poly(tractor_center, θ1, dp.tractor_length, dp.tractor_width)
    plot!(plt, tx, ty; lw=1, fill=(true, 0.10), label=false)

    # Trailer body: center around trailer axle forward by half trailer length
    trailer_center = P_tr + 0.5*dp.trailer_length * SVector(cos(θ2), sin(θ2))
    rx, ry = rect_poly(trailer_center, θ2, dp.trailer_length, dp.trailer_width)
    plot!(plt, rx, ry; lw=1, fill=(true, 0.08), label=false)

    # --- Axles (lines) ---
    if show_axes
        # rear axle line (perpendicular to θ1)
        n1 = SVector(-sin(θ1), cos(θ1))
        A1a = P_rear + dp.axle_halfwidth*n1
        A1b = P_rear - dp.axle_halfwidth*n1
        plot!(plt, [A1a[1], A1b[1]], [A1a[2], A1b[2]]; lw=2, label=false)

        # front axle line
        A2a = P_front + dp.axle_halfwidth*n1
        A2b = P_front - dp.axle_halfwidth*n1
        plot!(plt, [A2a[1], A2b[1]], [A2a[2], A2b[2]]; lw=2, label=false)

        # trailer axle line (perp to θ2)
        n2 = SVector(-sin(θ2), cos(θ2))
        Ta = P_tr + dp.axle_halfwidth*n2
        Tb = P_tr - dp.axle_halfwidth*n2
        plot!(plt, [Ta[1], Tb[1]], [Ta[2], Tb[2]]; lw=2, label=false)
    end

    # --- Wheels (rectangles) ---
    # rear wheels aligned with θ1
    n1 = SVector(-sin(θ1), cos(θ1))
    for s in (-1.0, 1.0)
        c = P_rear + s*dp.axle_halfwidth*n1
        wx, wy = rect_poly(c, θ1, dp.wheel_length, dp.wheel_width)
        plot!(plt, wx, wy; lw=1, fill=(true, 0.25), label=false)
    end

    # front wheels aligned with θ1+δ
    θw = θ1 + δ
    n1f = SVector(-sin(θ1), cos(θ1))
    for s in (-1.0, 1.0)
        c = P_front + s*dp.axle_halfwidth*n1f
        wx, wy = rect_poly(c, θw, dp.wheel_length, dp.wheel_width)
        plot!(plt, wx, wy; lw=1, fill=(true, 0.25), label=false)
    end

    # trailer wheels aligned with θ2
    n2 = SVector(-sin(θ2), cos(θ2))
    for s in (-1.0, 1.0)
        c = P_tr + s*dp.axle_halfwidth*n2
        wx, wy = rect_poly(c, θ2, dp.wheel_length, dp.wheel_width)
        plot!(plt, wx, wy; lw=1, fill=(true, 0.25), label=false)
    end

    # --- Hitch link ---
    plot!(plt, [P_rear[1], H[1]], [P_rear[2], H[2]]; lw=2, label=false)

    # --- Heading arrows ---
    if show_heading
        a1 = P_rear + 1.2*SVector(cos(θ1), sin(θ1))
        plot!(plt, [P_rear[1], a1[1]], [P_rear[2], a1[2]]; lw=2, label=false)
        a2 = P_tr + 1.2*SVector(cos(θ2), sin(θ2))
        plot!(plt, [P_tr[1], a2[1]], [P_tr[2], a2[2]]; lw=2, label=false)

        # steering direction indicator at front axle
        sdir = P_front + 1.0*SVector(cos(θw), sin(θw))
        plot!(plt, [P_front[1], sdir[1]], [P_front[2], sdir[2]]; lw=2, ls=:dash, label=false)
    end

    # --- φ arc (simple polyline arc around hitch) ---
    if show_phi_arc
        ϕ = x[4]
        r = 1.0
        t0 = θ1
        ts = range(t0, t0 + ϕ; length=25)
        ax = [H[1] + r*cos(t) for t in ts]
        ay = [H[2] + r*sin(t) for t in ts]
        plot!(plt, ax, ay; lw=2, ls=:dot, label=false)
    end

    return plt
end

function plot_xy_obstacles!(plt, obs2d; alpha=0.25)
    for ob in obs2d
        x1l,y1l = ob.lb[1], ob.lb[2]
        x1u,y1u = ob.ub[1], ob.ub[2]
        xs = [x1l,x1u,x1u,x1l,x1l]
        ys = [y1l,y1l,y1u,y1u,y1l]
        plot!(plt, xs, ys; lw=1, fill=(true, alpha), label=false)
    end
    return plt
end

function live_vehicle_progression(
    p, dp, traj::ST.Control_trajectory, xl, yl;
    domain=nothing,
    obstacles2d=nothing,
    every=1,
    dt=0.05,
    giffile::Union{Nothing,String}=nothing,
    fps::Int=20,
)
    states = traj.states.seq
    inputs = traj.inputs.seq

    xs = [x[1] for x in states]
    ys = [x[2] for x in states]

    # --- GIF MODE ---
    if giffile !== nothing
        anim = @animate for k in 1:every:length(states)
            plt = plot(; aspect_ratio=:equal, xlims=xl, ylims=yl, legend=false, size=(700,700))
            if domain !== nothing
                plot!(plt, domain; color=:grey, opacity=0.1)
            end
            if obstacles2d !== nothing
                plot_xy_obstacles!(plt, obstacles2d, color=:black)
            end
            plot!(plt, xs, ys; lw=1)
            uk = (k <= length(inputs)) ? inputs[k] : inputs[end]
            draw_articulated!(plt, p, dp, states[k], uk)
        end

        gif(anim, giffile; fps=fps)
        return anim
    end

    # --- LIVE MODE ---
    for k in 1:every:length(states)
        plt = plot(; aspect_ratio=:equal, xlims=xl, ylims=yl, legend=false, size=(700,700))
        if domain !== nothing
            plot!(plt, domain; color=:grey, opacity=0.1)
        end
        if obstacles2d !== nothing
            plot_xy_obstacles!(plt, obstacles2d, color=:black)
        end
        plot!(plt, xs, ys; lw=1)
        uk = (k <= length(inputs)) ? inputs[k] : inputs[end]
        draw_articulated!(plt, p, dp, states[k], uk)
        display(plt)
        sleep(dt)
    end

    return nothing
end


end # module
