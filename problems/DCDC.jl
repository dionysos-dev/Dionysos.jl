module DCDC
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
const CO = DI.Control
const SY = DI.Symbolic

function dynamicofsystem(vs = 1.0, rL = 0.05, xL = 3.0, rC = 0.005, xC = 70.0, r0 = 1.0, ngrowthbound = 5)
    # Definition of the dynamics functions $f_p$ of the system:
    b = SVector(vs/xL, 0.0);
    A1 = SMatrix{2,2}(-rL/xL, 0.0, 0.0, -1.0/xC/(r0+rC))
    A2 = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
        -r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC))
    F_sys = let b = b, A1 = A1, A2 = A2
        (x, u) -> u[1] == 1 ? A1*x + b : A2*x + b
    end
    A2_abs = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
                        r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC))
    L_growthbound = let A1 = A1, A2_abs = A2_abs
        u -> u[1] == 1 ? A1 : A2_abs
    end
    return F_sys, L_growthbound, ngrowthbound
end


function system(
    sysnoise = SVector(0.0, 0.0),
    measnoise = SVector(0.0, 0.0),
    tstep = 0.5,
    nsys = 5,
    _X_ = UT.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85)),
    _U_ = UT.HyperRectangle(SVector(1), SVector(2)),
    xdim=2,
    udim=1
)
    # Definition of the dynamics functions $f_p$ of the system:
    F_sys,L_growthbound,ngrowthbound=dynamicofsystem();
    contsys = ST.NewControlSystemGrowthRK4(tstep, F_sys, L_growthbound, sysnoise,
                                        measnoise, nsys, ngrowthbound);
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(contsys, xdim, udim, _X_, _U_)
end

function problem()
    sys = system()
    return PB.SafetyProblem(sys, sys.X, sys.X, PB.Infinity())
end

end
