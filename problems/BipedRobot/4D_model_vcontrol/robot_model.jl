module RobotModel

using StaticArrays
using Plots
using Dionysos
using MathematicalSystems

const UT = Dionysos.Utils
const PB = Dionysos.Problem
const ST = Dionysos.System

"""
States (angles) and inputs (angular velocities):
    1 -> left hip
    2 -> left knee
    3 -> right hip
    4 -> right knee
"""

robot_dynamics() = (x,u) -> SVector{4}(u[1], u[2], u[3], u[4])
robot_dynamics_dt(dt) = (x,u) -> x + dt*u

# For ẋ = u, ∂f/∂x = 0, so a constant zero matrix bound is valid.
# Returned function signature matches Dionysos' growth-bound constructors: L(u).
function jacobian_bound()
    return (u) -> @SMatrix [
        0.0 0.0;
        0.0 0.0
    ]
end

function system(;
    # theta_max = 30°
    # omega_max = 15°/s
    _X_ = UT.HyperRectangle(SVector(-pi/3, -pi/3, -pi/3, -pi/3), SVector(pi/3, pi/3, pi/3, pi/3)),
    _U_ = UT.HyperRectangle(SVector(-pi/6, -pi/6, -pi/6, -pi/6), SVector(pi/6, pi/6, pi/6, pi/6)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        robot_dynamics(),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function dt_system(
    dt;
    # theta_max = 30°
    # omega_max = 15°/s
    _X_ = UT.HyperRectangle(SVector(-pi/3, -pi/3, -pi/3, -pi/3), SVector(pi/3, pi/3, pi/3, pi/3)),
    _U_ = UT.HyperRectangle(SVector(-pi/6, -pi/6, -pi/6, -pi/6), SVector(pi/6, pi/6, pi/6, pi/6))
)
    return MathematicalSystems.ConstrainedBlackBoxControlDiscreteSystem(
        robot_dynamics_dt(dt),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function optimal_control_problem(;
    _X_ = UT.HyperRectangle(SVector(-pi/3, -pi/3, -pi/3, -pi/3), SVector(pi/3, pi/3, pi/3, pi/3)),
    _U_ = UT.HyperRectangle(SVector(-pi/6, -pi/6, -pi/6, -pi/6), SVector(pi/6, pi/6, pi/6, pi/6)),
    _I_ = UT.HyperRectangle(SVector(-pi/12, -pi/12, -pi/12, -pi/12), SVector(pi/12, pi/12, pi/12, pi/12)),
    _T_ = UT.HyperRectangle(SVector(-pi/3+pi/12, -pi/12, pi/3-pi/12, -pi/12), SVector(-pi/3, pi/12, pi/3, pi/12)),
)
    sys = system(; _X_ = _X_, _U_ = _U_)
    return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
end

end