module SimplePendulum

using StaticArrays
using MathematicalSystems
using Dionysos
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const PB = Dionysos.Problem
const ST = Dionysos.System

function dynamic(;l = 1.0, g = 9.81)
    return (x, u) -> begin
        return SVector{2}(x[2], -(g / l) * sin(x[1]) + u[1])
    end
end

function jacobian_bound(l = 1.0, g = 9.81)
    function L_growthbound(u)
        return SMatrix{2, 2}(0.0, 1.0, (g / l), 0)
    end
end

function system(;
    l = 1.0, g = 9.81,
    _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0)),
    _U_ = UT.HyperRectangle(SVector(-6), SVector(6)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(l=l, g=g),
        Dionysos.Utils.get_dims(_X_),
        Dionysos.Utils.get_dims(_U_),
        _X_,
        _U_,
    )
end

function optimal_control_problem(; 
        l = 1.0, g = 9.81,
        _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0)),
        _U_ = UT.HyperRectangle(SVector(-6), SVector(6)),
        _I_ = UT.HyperRectangle(SVector(-5.0 * pi / 180.0, -0.5),SVector(5.0 * pi / 180.0, 0.5)),
        _T_ = UT.HyperRectangle(SVector(pi - 10.0 * pi / 180.0, -3.0), SVector(pi + 10.0 * pi / 180.0, 3.0))
        )
        sys = system(l = l, g = g, _X_ = _X_, _U_ = _U_)
        return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
end
function safety_control_problem(; 
        l = 1.0, g = 9.81,
        _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0)),
        _U_ = UT.HyperRectangle(SVector(-6), SVector(6)),
        _I_ = UT.HyperRectangle(SVector(pi - 5.0 * pi / 180.0, -0.5), SVector(pi + 5.0 * pi / 180.0, 0.5)),
        _S_ = UT.HyperRectangle(SVector(pi - 8.0 * pi / 180.0, -1.0), SVector(pi + 8.0 * pi / 180.0, 1.0))
    )
        sys = system(l = l, g = g, _X_ = _X_, _U_ = _U_)
        return PB.SafetyControlProblem(sys, _I_, _S_, PB.Infinity())
end

end
