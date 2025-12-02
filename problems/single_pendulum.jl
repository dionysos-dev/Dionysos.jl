module SinglePendulum

using StaticArrays
using MathematicalSystems, HybridSystems
import Dionysos
using Dionysos.Problem
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

g = 9.81
l = 1.0

function dynamic()
    return (x::SVector{2, Float64}, u::SVector{1, Float64}) -> begin
        return SVector{2}(
            x[2], 
            u[1] - g / l * sin(x[1]),
            )
    end
end

function jacobian_bound(u)
    return SMatrix{2, 2}(
        0.0, 1.0, 
        (g / l), 0,
        )
end

function system(_X_, _U_)

    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(),
        Dionysos.Utils.get_dims(_X_),
        Dionysos.Utils.get_dims(_U_),
        _X_,
        _U_,
    )
end

function problem(; _U_ = UT.HyperRectangle(SVector(-3.0), SVector(3.0)), transition_cost = UT.ConstantControlFunction(1.0))

    _X_ = UT.HyperRectangle(SVector(-pi, -10.0), SVector(2pi, 10.0))
    _I_ = UT.HyperRectangle(SVector(-5.0 * pi / 180.0, -0.5), SVector(5.0 * pi / 180.0, 0.5))
    _T_ = UT.HyperRectangle(SVector(pi - 30.0 * pi / 180.0, -1.0), SVector(pi + 30.0 * pi / 180.0, 1.0))
    sys = system(_X_, _U_)
    problem = OptimalControlProblem(
        sys,
        _I_,
        _T_,
        UT.ZeroFunction(),
        transition_cost,
        PR.Infinity(),
    )
end 

end # Module SinglePendulum