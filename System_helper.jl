# Dynamics F 
function _F_()
    return (x, u) -> begin
        α = atan(tan(u[2]) / 2)
        return SVector{3}(
            u[1] * cos(α + x[3]) / cos(α),
            u[1] * sin(α + x[3]) / cos(α),
            u[1] * tan(u[2]),
        )
    end
end


function jacobian_bound()
    return u -> begin
        β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
        return SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
    end
end



function get_obstacles(_X_, X1_lb, X1_ub, X2_lb, X2_ub; margin::Float64 = 0.0)
    return [
        UT.HyperRectangle(
            SVector(
                max(_X_.lb[1], x1lb - margin),
                max(_X_.lb[2], x2lb - margin),
                _X_.lb[3],
            ),
            SVector(
                min(_X_.ub[1], x1ub + margin),
                min(_X_.ub[2], x2ub + margin),
                _X_.ub[3],
            ),
        )
        for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
    ]
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



function system(_X_; _U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0)))
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        _F_(),
        Dionysos.Utils.get_dims(_X_),
        Dionysos.Utils.get_dims(_U_),
        _X_,
        _U_,
    )
end



function problem(x1_lb, x1_ub, x2_lb, x2_ub; simple = false, transition_cost = nothing)
    # 1) Base box X_box
    if simple
        X_box = UT.HyperRectangle(
            SVector(0.0, 0.0, -pi - 0.4),
            SVector(10.0, 10.0, pi + 0.4),
        )
    else
        X_box = UT.HyperRectangle(
            SVector(0.0, 0.0, -pi - 0.4),
            SVector(10.0, 10.0, pi + 0.4),
        )
    end
    U_box = U_BOX

    # 2) Obstacles in that box
    obs = get_obstacles(X_box, x1_lb, x1_ub, x2_lb, x2_ub; margin = OBSTACLE_MARGIN)
    obstacles_LU = UT.LazyUnionSetArray(obs)

    # 3) Free state space for the system
    X_free = UT.LazySetMinus(X_box, obstacles_LU)

    # 4) System uses X_free, but EmptyProblem stores X_box
    sys = system(X_free; _U_ = U_box)

    return PB.EmptyProblem(sys, X_box)
end


x1_lb = [0.5, 0.5, 3.5, 3.5, 4.0, 5.5]
x1_ub = [2.5, 2.5, 4.0, 4.0, 6.0, 6.0]
x2_lb = [5.0, 6.5, 1.0, 4.0, 4.0, 4.5]
x2_ub = [5.5, 9.0, 3.0, 6.0, 4.5, 6.0]
