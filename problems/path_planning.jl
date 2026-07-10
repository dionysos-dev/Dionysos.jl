module PathPlanning

using StaticArrays
using MathematicalSystems, HybridSystems
using Plots
import LazySets

using Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem

# System eq x' = F_sys(x, u)
function dynamic()
    return (x, u) -> begin
        α = atan(tan(u[2]) / 2)
        return SVector{3}(
            u[1] * cos(α + x[3]) / cos(α),
            u[1] * sin(α + x[3]) / cos(α),
            u[1] * tan(u[2]),
        )
    end
end

# We define the growth bound function of $f$:
function jacobian_bound()
    return u -> begin
        β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
        return SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
    end
end

function jacobian()
    return (x, u) -> begin
        α = atan(tan(u[2]) / 2)
        β = u[1] / cos(α)
        return SMatrix{3, 3}(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -β * sin(α + x[3]),
            β * cos(α + x[3]),
            0.0,
        )
    end
end

function bound_norm_jacobian()
    return u -> abs(u[1] / cos(atan(tan(u[2]) / 2)))
end

function bound_norm_hessian_tensor()
    return u -> abs(u[1] / cos(atan(tan(u[2]) / 2)))
end

function filter_obstacles(_X_, _I_, _T_, obs)
    obstacles = typeof(_X_)[]
    for ob in obs
        if ob ⊆ _X_ && isempty(ob ∩ _I_) && isempty(ob ∩ _T_)
            push!(obstacles, ob)
        end
    end
    obstacles_LU = UT.set_union(obstacles)
    return obstacles_LU
end

function get_obstacles(
    _X_;
    X1_lb = [1.0, 2.2, 2.2, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3],
    X1_ub = [1.2, 2.4, 2.4, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0],
    X2_lb = [0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2],
    X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6],
)
    return [
        UT.box(
            SVector(x1lb, x2lb, LazySets.low(_X_, 3)),
            SVector(x1ub, x2ub, LazySets.high(_X_, 3)),
        ) for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
    ]
end

function system(_X_; _U_ = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0)))
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(),
        Dionysos.Utils.get_dim(_X_),
        Dionysos.Utils.get_dim(_U_),
        _X_,
        _U_,
    )
end

""""
    problem()

This function create the system with `PathPlanning.system`.

Then, we define initial and target domains for the state of the system.

Finally, we instantiate our Reachability Problem as an OptimalControlProblem 
with the system, the initial and target domains, and null cost functions.
"""
function problem(; simple = false, transition_cost = nothing)
    if simple
        _X_ = UT.box(SVector(0.0, 0.0, -pi - 0.4), SVector(4.0, 10.0, pi + 0.4))
        _I_ = UT.box(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
        _T_ = UT.box(SVector(3.0, 0.3, -100.0), SVector(3.6, 0.8, 100.0))
    else
        _X_ = UT.box(SVector(-0.1, -0.1, -pi - 0.4), SVector(10.1, 10.1, pi + 0.4))
        _I_ = UT.box(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
        _T_ = UT.box(SVector(9.0, 0.3, -100.0), SVector(9.6, 0.8, 100.0))
    end
    obs = get_obstacles(_X_)
    obstacles_LU = filter_obstacles(_X_, _I_, _T_, obs)
    _X_ = UT.set_minus(_X_, obstacles_LU)
    sys = system(_X_)

    problem = PR.OptimalControlProblem(sys, _I_, _T_, nothing, transition_cost)
    return problem
end

function plot_xy_obstacles!(fig, obstacles; alpha = 0.25)
    for ob in obstacles
        x1l, x2l = LazySets.low(ob, 1), LazySets.low(ob, 2)
        x1u, x2u = LazySets.high(ob, 1), LazySets.high(ob, 2)

        xs = [x1l, x1u, x1u, x1l, x1l]
        ys = [x2l, x2l, x2u, x2u, x2l]

        plot!(fig, xs, ys; linewidth = 1, fill = (true, alpha), color = :black, label = "")
    end

    return fig
end

function system_plot!(;
    car_length = 0.45,
    car_width = 0.25,
    arrow_length = 0.6,
    xlims = (-0.1, 10.1),
    ylims = (-0.1, 10.1),
    obstacles = nothing,
    show_input = true,
)
    return function (fig, x, u)
        px = Float64(x[1])
        py = Float64(x[2])
        θ = Float64(x[3])

        v = Float64(u[1])
        δ = Float64(u[2])

        if obstacles !== nothing
            plot_xy_obstacles!(fig, obstacles)
        end

        R = @SMatrix [
            cos(θ) -sin(θ)
            sin(θ) cos(θ)
        ]

        c = SVector(px, py)

        corners_body = [
            SVector(car_length / 2, car_width / 2),
            SVector(car_length / 2, -car_width / 2),
            SVector(-car_length / 2, -car_width / 2),
            SVector(-car_length / 2, car_width / 2),
            SVector(car_length / 2, car_width / 2),
        ]

        body = [c + R * q for q in corners_body]

        plot!(
            fig,
            [p[1] for p in body],
            [p[2] for p in body];
            linewidth = 2,
            fill = (true, 0.25),
            color = :blue,
            label = "",
        )

        # heading arrow
        p_head = c + arrow_length * SVector(cos(θ), sin(θ))

        plot!(
            fig,
            [px, p_head[1]],
            [py, p_head[2]];
            linewidth = 3,
            color = :blue,
            label = "",
        )

        scatter!(fig, [px], [py]; markersize = 5, color = :blue, label = "")

        if show_input
            annotate!(
                fig,
                xlims[1] + 0.05 * (xlims[2] - xlims[1]),
                ylims[1] + 0.05 * (ylims[2] - ylims[1]),
                text("v = $(round(v; digits = 2))\nδ = $(round(δ; digits = 2))", 10),
            )
        end

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end
