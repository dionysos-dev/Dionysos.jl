module Pendulum
using Test

using StaticArrays
using MathematicalSystems

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic

function dynamicofsystem(l = 1.0, g = 9.81)
    function F_sys(x, u)
        return SVector{2}(x[2], -(g / l) * sin(x[1]) + u[1])
    end

    function L_growthbound(u)
        return SMatrix{2, 2}(0.0, 1.0, (g / l), 0)
    end

    return F_sys, L_growthbound
end

function system(;
    sysnoise = SVector(0.0, 0.0),
    measnoise = SVector(0.0, 0.0),
    tstep = 0.1,
    nsys = 5,
    _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π + pi, 5.0)),
    _U_ = UT.HyperRectangle(SVector(-6), SVector(6)), #-5 -11
    xdim = 2,
    udim = 1,
)
    F_sys, L_growthbound = dynamicofsystem()
    ngrowthbound = 5
    contsys = ST.NewControlSystemGrowthRK4(
        tstep,
        F_sys,
        L_growthbound,
        sysnoise,
        measnoise,
        nsys,
        ngrowthbound,
    )
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        contsys,
        xdim,
        udim,
        _X_,
        _U_,
    )
end

function problem(; objective = "reachability-up-low_power")
    if objective == "safety-up"
        sys = system()
        _I_ = UT.HyperRectangle(
            SVector(pi - 5.0 * pi / 180.0, -0.5),
            SVector(pi + 5.0 * pi / 180.0, 0.5),
        )
        _S_ = UT.HyperRectangle(
            SVector(pi - 8.0 * pi / 180.0, -1.0),
            SVector(pi + 8.0 * pi / 180.0, 1.0),
        )
        return PB.SafetyProblem(sys, _I_, _S_, PB.Infinity())
    elseif objective == "safety-down"
        sys = system()
        _I_ = UT.HyperRectangle(
            SVector(-5.0 * pi / 180.0, -0.5),
            SVector(5.0 * pi / 180.0, 0.5),
        )
        _S_ = UT.HyperRectangle(
            SVector(-8.0 * pi / 180.0, -1.0),
            SVector(8.0 * pi / 180.0, 1.0),
        )
        return PB.SafetyProblem(sys, _I_, _S_, PB.Infinity())
    elseif objective == "reachability-up-ultra_low_power"
        sys = system(;
            _X_ = UT.HyperRectangle(SVector(-π, -10.0), SVector(π + pi, 10.0)),
            _U_ = UT.HyperRectangle(SVector(-3.0), SVector(3.0)),
        )
        _I_ = UT.HyperRectangle(
            SVector(-5.0 * pi / 180.0, -0.5),
            SVector(5.0 * pi / 180.0, 0.5),
        )
        _T_ = UT.HyperRectangle(
            SVector(pi - 5.0 * pi / 180.0, -1.0),
            SVector(pi + 5.0 * pi / 180.0, 1.0),
        )
        return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
    elseif objective == "reachability-up-low_power"
        sys = system(; _U_ = UT.HyperRectangle(SVector(-4), SVector(4)))
        _I_ = UT.HyperRectangle(
            SVector(-5.0 * pi / 180.0, -0.5),
            SVector(5.0 * pi / 180.0, 0.5),
        )
        _T_ = UT.HyperRectangle(
            SVector(pi - 5.0 * pi / 180.0, -1.0),
            SVector(pi + 5.0 * pi / 180.0, 1.0),
        )
        return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
    elseif objective == "reachability-up-high_power"
        sys = system(; _U_ = UT.HyperRectangle(SVector(-10), SVector(10)))
        _I_ = UT.HyperRectangle(
            SVector(-5.0 * pi / 180.0, -0.5),
            SVector(5.0 * pi / 180.0, 0.5),
        )
        _T_ = UT.HyperRectangle(
            SVector(pi - 5.0 * pi / 180.0, -2.0),
            SVector(pi + 5.0 * pi / 180.0, 2.0),
        )
        return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
    end
end

end
