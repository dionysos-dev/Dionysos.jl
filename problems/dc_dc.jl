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
const SY = DI.Symbolic

function A1(; xL = 3.0, xC = 70.0, r0 = 1.0, rL = 0.05, rC = 0.005)
    return SMatrix{2, 2}(-rL / xL, 0.0, 0.0, -1.0 / xC / (r0 + rC))
end

function A2(; xL = 3.0, xC = 70.0, r0 = 1.0, rL = 0.05, rC = 0.005)
    return SMatrix{2, 2}(
        -(rL + r0 * rC / (r0 + rC)) / xL,
        5.0 * r0 / (r0 + rC) / xC,
        -r0 / (r0 + rC) / xL / 5.0,
        -1.0 / xC / (r0 + rC),
    )
end

function A2_abs(; xL = 3.0, xC = 70.0, r0 = 1.0, rL = 0.05, rC = 0.005)
    return SMatrix{2, 2}(
        -(rL + r0 * rC / (r0 + rC)) / xL,
        5.0 * r0 / (r0 + rC) / xC,
        r0 / (r0 + rC) / xL / 5.0,
        -1.0 / xC / (r0 + rC),
    )
end

function jacobian_bound(; kws...)
    return let A1 = A1(; kws...), A2_abs = A2_abs(; kws...)
        u -> u[1] == 1 ? A1 : A2_abs
    end
end

function DF_sys(; kws...)
    return let A1 = A1(; kws...), A2 = A2(; kws...)
        u -> u[1] == 1 ? A1 : A2
    end
end

function dynamic(; vs = 1.0, xL = 3.0, kws...)
    # Definition of the dynamics functions $f_p$ of the system:
    b = SVector(vs / xL, 0.0)
    return let b = b, A1 = A1(; xL, kws...), A2 = A2(; xL, kws...)
        (x, u) -> u[1] == 1 ? A1 * x + b : A2 * x + b
    end
end

function system(;
    _X_ = UT.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85)),
    _U_ = UT.HyperRectangle(SVector(1), SVector(2)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(),
        Dionysos.Utils.get_dims(_X_),
        Dionysos.Utils.get_dims(_U_),
        _X_,
        _U_,
    )
end

function problem()
    sys = system()
    _I_ = UT.HyperRectangle(SVector(1.19, 5.59), SVector(1.21, 5.61))
    return PB.SafetyProblem(sys, _I_, sys.X, PB.Infinity())
end

end
