module ToyProblem

using StaticArrays
using MathematicalSystems
using Dionysos

const UT = Dionysos.Utils
const PB = Dionysos.Problem

# 2D integrator: ẋ = u
dynamic() = (x, u) -> SVector{2}(u[1], u[2])

# For ẋ = u, ∂f/∂x = 0, so a constant zero matrix bound is valid.
# Returned function signature matches Dionysos' growth-bound constructors: L(u).
function jacobian_bound()
    return (u) -> @SMatrix [
        0.0 0.0;
        0.0 0.0
    ]
end

function system(;
    _X_ = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0)),
    _U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function optimal_control_problem(;
    _X_ = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0)),
    _U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0)),
    _I_ = UT.HyperRectangle(SVector(-1.6, -1.6), SVector(-1.4, -1.4)),
    _T_ = UT.HyperRectangle(SVector(-0.2, -0.2), SVector(0.2, 0.2)),
)
    sys = system(; _X_ = _X_, _U_ = _U_)
    return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
end

end # module
