module SimplePendulum

using StaticArrays
using MathematicalSystems
using Dionysos
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const PB = Dionysos.Problem
const ST = Dionysos.System

function dynamic(; l = 1.0, g = 9.81)
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
    l = 1.0,
    g = 9.81,
    _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0)),
    _U_ = UT.HyperRectangle(SVector(-6), SVector(6)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(; l = l, g = g),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function safety_problem(;
    l = 1.0,
    g = 9.81,
    _X_ = UT.HyperRectangle(SVector(-π, -1.5), SVector(π, 1.5)),
    _U_ = UT.HyperRectangle(SVector(-4.0), SVector(4.0)),
    _I_ = UT.HyperRectangle(SVector(pi-3.0*pi/180.0, -0.5), SVector(pi+3.0*pi/180.0, 0.5)),
    _S_ = UT.HyperRectangle(
        SVector(pi-15.0*pi/180.0, -1.0),
        SVector(pi+15.0*pi/180.0, 1.0),
    ),
    objective = "safety_up",
)
    if objective == "safety_down"
        _X_ = UT.HyperRectangle(SVector(-π, -1.5), SVector(π, 1.5))
        _U_ = UT.HyperRectangle(SVector(-4.0), SVector(4.0))
        _I_ = UT.HyperRectangle(SVector(-3.0*pi/180.0, -0.5), SVector(3.0*pi/180.0, 0.5))
        _S_ = UT.HyperRectangle(SVector(-15.0*pi/180.0, -1.0), SVector(15.0*pi/180.0, 1.0))
    elseif objective == "safety_up"
        _X_ = UT.HyperRectangle(SVector(-π, -1.5), SVector(π, 1.5))
        _U_ = UT.HyperRectangle(SVector(-4.0), SVector(4.0))
        _I_ =
            UT.HyperRectangle(SVector(pi-3.0*pi/180.0, -0.5), SVector(pi+3.0*pi/180.0, 0.5))
        _S_ = UT.HyperRectangle(
            SVector(pi-15.0*pi/180.0, -1.0),
            SVector(pi+15.0*pi/180.0, 1.0),
        )
    end
    sys = system(; l = l, g = g, _X_ = _X_, _U_ = _U_)
    return PB.SafetyProblem(sys, _I_, _S_, PB.Infinity())
end

function optimal_control_problem(;
    l = 1.0,
    g = 9.81,
    _X_ = UT.HyperRectangle(SVector(-π, -7.0), SVector(π, 7.0)),
    _U_ = UT.LazySetMinus(
        UT.HyperRectangle(SVector(-2.5), SVector(2.5)),
        UT.HyperRectangle(SVector(-0.5), SVector(0.5)),
    ),
    _I_ = UT.HyperRectangle(
        SVector(-5.0 * pi / 180.0, -0.2),
        SVector(5.0 * pi / 180.0, 0.2),
    ),
    _T_ = UT.HyperRectangle(
        SVector(pi - 15.0 * pi / 180.0, -1.0),
        SVector(pi + 15.0 * pi / 180.0, 1.0),
    ),
    _O_ = UT.HyperRectangle(
        SVector(-pi + 16.0 * pi / 180.0, -7.0),
        SVector(-pi + 38.0 * pi / 180.0, 7.0),
    ),
    objective = "reachability_up_low_power",
)
    if objective == "reachability_up_high_power"
        _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0))
        _U_ = UT.LazySetMinus(
            UT.HyperRectangle(SVector(-10.0), SVector(10.0)),
            UT.HyperRectangle(SVector(-0.5), SVector(0.5)),
        )
        _I_ = UT.HyperRectangle(
            SVector(-5.0 * pi / 180.0, -0.2),
            SVector(5.0 * pi / 180.0, 0.2),
        )
        _T_ = UT.HyperRectangle(
            SVector(pi - 15.0 * pi / 180.0, -1.0),
            SVector(pi + 15.0 * pi / 180.0, 1.0),
        )
    elseif objective == "reachability_up_medium_power"
        _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π, 5.0))
        _U_ = UT.LazySetMinus(
            UT.HyperRectangle(SVector(-7.0), SVector(7.0)),
            UT.HyperRectangle(SVector(-0.5), SVector(0.5)),
        )
        _I_ = UT.HyperRectangle(
            SVector(-5.0 * pi / 180.0, -0.2),
            SVector(5.0 * pi / 180.0, 0.2),
        )
        _T_ = UT.HyperRectangle(
            SVector(pi - 15.0 * pi / 180.0, -1.0),
            SVector(pi + 15.0 * pi / 180.0, 1.0),
        )
    elseif objective == "reachability_up_low_power"
        _X_ = UT.HyperRectangle(SVector(-π, -7.0), SVector(π, 7.0))
        _U_ = UT.LazySetMinus(
            UT.HyperRectangle(SVector(-2.5), SVector(2.5)),
            UT.HyperRectangle(SVector(-0.5), SVector(0.5)),
        )
        _I_ = UT.HyperRectangle(
            SVector(-5.0 * pi / 180.0, -0.2),
            SVector(5.0 * pi / 180.0, 0.2),
        )
        _T_ = UT.HyperRectangle(
            SVector(pi - 15.0 * pi / 180.0, -1.0),
            SVector(pi + 15.0 * pi / 180.0, 1.0),
        )
    end
    _X_ = _O_ !== nothing ? UT.LazySetMinus(_X_, _O_) : _X_
    sys = system(; l = l, g = g, _X_ = _X_, _U_ = _U_)
    return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
end

end
