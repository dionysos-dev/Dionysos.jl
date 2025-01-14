module Pendulum
using Test
# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).
using StaticArrays
using MathematicalSystems
# At this point, we import the useful Dionysos sub-module for this problem:
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic
function dynamicofsystem(
    l = 1.0,
    g = 9.81,
)
    function F_sys(x, u)
        return SVector{2}(
            x[2],
            -(g/l)*sin(x[1]) + u[1],  # u = N / kg 
        )
    end
    function L_growthbound(u)
        return SMatrix{2, 2}(0.0, 1.0, (g/l), 0)
    end
    return F_sys, L_growthbound
end
function system(;
    sysnoise = SVector(0.0, 0.0),
    measnoise = SVector(0.0, 0.0),
    tstep = 0.1, #0.5
    nsys = 5,
    _X_ = UT.HyperRectangle(SVector(-π, -5.0), SVector(π+pi, 5.0)), 
    _U_ = UT.HyperRectangle(SVector(-6.0), SVector(6.0)),#11 #8 #6
    xdim = 2,
    udim = 1,
    approx_mode = "growth",
)
    # Definition of the dynamics functions $f_p$ of the system:
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
function problem(; approx_mode = "growth")
    sys = system(; approx_mode = approx_mode)
    ## stable equilibrium
    _I_ = UT.HyperRectangle(SVector(-5.0*pi/180.0, -0.5), SVector(5.0*pi/180.0, 0.5))
    _S_ = UT.HyperRectangle(SVector(-30.0*pi/180.0, -1.0), SVector(30.0*pi/180.0, 1.0))
    _T_ = UT.HyperRectangle(SVector(pi-30.0*pi/180.0, -1.0), SVector(pi+30.0*pi/180.0, 1.0))
    ## unstable equilibrium
    return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
    #return PB.SafetyProblem(sys, _I_, _S_, PB.Infinity())
end
end