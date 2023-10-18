module PathPlanning

using StaticArrays
using MathematicalSystems, HybridSystems
using Dionysos
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const PB = Dionysos.Problem
const ST = Dionysos.System

function dynamicofsystem()
    # System eq x' = F_sys(x, u)
    function F_sys(x, u)
        α = atan(tan(u[2]) / 2)
        return SVector{3}(
            u[1] * cos(α + x[3]) / cos(α),
            u[1] * sin(α + x[3]) / cos(α),
            u[1] * tan(u[2]),
        )
    end

    # We define the growth bound function of $f$:
    ngrowthbound = 5
    function L_growthbound(u)
        β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
        return SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
    end
    function DF_sys(x, u)
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
    bound_DF(u) = abs(u[1] / cos(atan(tan(u[2]) / 2)))
    bound_DDF(u) = abs(u[1] / cos(atan(tan(u[2]) / 2)))

    return F_sys, L_growthbound, ngrowthbound, DF_sys, bound_DF, bound_DDF
end

function _initTargetSets()
    _I_ = UT.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
    _T_ = UT.HyperRectangle(SVector(3.0, 0.3, -100.0), SVector(3.6, 0.8, 100.0))
    return _I_, _T_
end

function filter_obstacles(_X_, obs)
    _I_, _T_ = _initTargetSets()
    obstacles = typeof(_X_)[]
    for ob in obs
        if ob ⊆ _X_ && isempty(ob ∩ _I_) && isempty(ob ∩ _T_)
            push!(obstacles, ob)
        end
    end
    obstacles_LU = UT.LazyUnionSetArray(obstacles)
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
        UT.HyperRectangle(SVector(x1lb, x2lb, _X_.lb[3]), SVector(x1ub, x2ub, _X_.ub[3]))
        for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
    ]
end

function system(
    _X_;
    _U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0)),
    xdim = 3,
    udim = 2,
    sysnoise = SVector(0.0, 0.0, 0.0),
    measnoise = SVector(0.0, 0.0, 0.0),
    tstep = 0.3,
    nsys = 5,
    approx_mode = "growth",
)
    F_sys, L_growthbound, ngrowthbound, DF_sys, bound_DF, bound_DDF = dynamicofsystem()
    contsys = nothing
    if approx_mode == "growth"
        contsys = ST.NewControlSystemGrowthRK4(
            tstep,
            F_sys,
            L_growthbound,
            sysnoise,
            measnoise,
            nsys,
            ngrowthbound,
        )
    elseif approx_mode == "linearized"
        contsys = ST.NewControlSystemLinearizedRK4(
            tstep,
            F_sys,
            DF_sys,
            bound_DF,
            bound_DDF,
            measnoise,
            nsys,
        )
    end
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        contsys,
        xdim,
        udim,
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
function problem(; simple = false, approx_mode = "growth")
    if simple
        _X_ = UT.HyperRectangle(SVector(0.0, 0.0, -pi - 0.4), SVector(4.0, 10.0, pi + 0.4))
        _I_ = UT.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
        _T_ = UT.HyperRectangle(SVector(3.0, 0.3, -100.0), SVector(3.6, 0.8, 100.0))
    else
        _X_ = UT.HyperRectangle(SVector(0.0, 0.0, -pi - 0.4), SVector(10.0, 10.0, pi + 0.4))
        _I_ = UT.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
        _T_ = UT.HyperRectangle(SVector(9.0, 0.3, -100.0), SVector(9.6, 0.8, 100.0))
    end
    obs = get_obstacles(_X_)
    obstacles_LU = filter_obstacles(_X_, obs)
    _X_ = UT.LazySetMinus(_X_, obstacles_LU)
    sys = system(_X_; approx_mode = approx_mode)
    problem = PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
    return problem
end

end
