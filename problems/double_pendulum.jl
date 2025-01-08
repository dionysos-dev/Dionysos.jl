module DoublePendulum
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
    l1 = 1.0,
    l2 = 1.0,
    m1=1.0,
    m2=1.0,
    g = 9.81,
)
    function F_sys(x, u)
        M = m1+m2
        Δθ = x[1]-x[2]
        α = m1+m2*sin(Δθ)^2
        return SVector{4}(
            x[3],
            x[4],
            -sin(Δθ)*(m2*l1*x[3]^2*cos(Δθ)+m2*l2*x[4]^2) - g*(M*sin(x[1])-m2*sin(x[2])*cos(Δθ))/(l1*α) + u[1],
            sin(Δθ)*(M*l1*x[3]^2+m2*l2*x[4]^2*cos(Δθ)) + g*(M*sin(x[1])*cos(Δθ)-M*sin(x[2]))/(l2*α),
        )
    end
    function L_growthbound(u)
        return 0.
    end
    return F_sys, L_growthbound
end
function system(;
    sysnoise = SVector(0.0, 0.0, 0.0, 0.0),
    measnoise = SVector(0.0, 0.0, 0.0, 0.0),
    tstep = 0.1, # time step du truc, comme max je prends 5 fois ça
    nsys = 5,
    _X_ = UT.HyperRectangle(SVector(-π, -π, -3, -3), SVector(π, π, 3, 3)),
    _U_ = UT.HyperRectangle(SVector(-11), SVector(11)),
    xdim = 4,
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
    _I_ = UT.HyperRectangle(SVector(-5.0*pi/180.0, -5.0*pi/180.0, -0.5, -0.5), SVector(5.0*pi/180.0, 5.0*pi/180.0, 0.5, 0.5))
    #_S_ = UT.HyperRectangle(SVector(-30.0*pi/180.0, -1.0), SVector(30.0*pi/180.0, 1.0))
    _S_ = UT.HyperRectangle(SVector(-30.0*pi/180.0, -30.0*pi/180.0, -1.0, -1.0), SVector(30.0*pi/180.0, 30.0*pi/180.0, 1.0, 1.0))

    return PB.SafetyProblem(sys, _I_, _S_, PB.Infinity())
end
end