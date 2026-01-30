module DoublePendulum

using StaticArrays
using MathematicalSystems
using Dionysos
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const PB = Dionysos.Problem
const ST = Dionysos.System

function dynamic(; l1 = 1.0, l2 = 1.0, m1 = 1.0, m2 = 1.0, g = 9.81)
    function f(x, u)
        M = m1 + m2
        Δθ = x[1] - x[2]
        α = m1 + m2 * sin(Δθ)^2
        return SVector{4}(
            x[3],
            x[4],
            -sin(Δθ) * (m2 * l1 * x[3]^2 * cos(Δθ) + m2 * l2 * x[4]^2) -
            g * (M * sin(x[1]) - m2 * sin(x[2]) * cos(Δθ)) / (l1 * α) + u[1],
            sin(Δθ) * (M * l1 * x[3]^2 + m2 * l2 * x[4]^2 * cos(Δθ)) +
            g * (M * sin(x[1]) * cos(Δθ) - M * sin(x[2])) / (l2 * α),
        )
    end
    return f
end

function jacobian_bound(l = 1.0, g = 9.81)
    function L_growthbound(u)
        return SMatrix{4, 4}(0.0, 1.0, (g / l), 0)
    end
end

function system(;
    l1 = 1.0,
    l2 = 1.0,
    m1 = 1.0,
    m2 = 1.0,
    g = 9.81,
    _X_ = UT.HyperRectangle(SVector(-π, -π, -5.0, -5.0), SVector(π, π, 5.0, 5.0)),
    _U_ = UT.HyperRectangle(SVector(-11.0), SVector(11.0)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(; l1 = l1, l2 = l2, m1 = m1, m2 = m2, g = g),
        UT.get_dims(_X_),
        UT.get_dims(_U_),
        _X_,
        _U_,
    )
end

function safety_problem(;
    l1 = 1.0,
    l2 = 1.0,
    m1 = 1.0,
    m2 = 1.0,
    g = 9.81,
    objective = "safety_down",
)
    if objective == "safety_down"
        _X_ = UT.HyperRectangle(
            SVector(-π/4.0, -π/4.0, -5.0, -5.0),
            SVector(π/4.0, π/4.0, 5.0, 5.0),
        )
        _U_ = UT.HyperRectangle(SVector(-2.5), SVector(2.5))
        _I_ = UT.HyperRectangle(
            SVector(-3.0*pi/180.0, -3.0*pi/180.0, -0.5, -0.5),
            SVector(3.0*pi/180.0, 3.0*pi/180.0, 0.5, 0.5),
        )
        _S_ = UT.HyperRectangle(
            SVector(-35.0*pi/180.0, -35.0*pi/180.0, -1.0, -1.0),
            SVector(35.0*pi/180.0, 35.0*pi/180.0, 1.0, 1.0),
        )
    end
    sys = system(; l1 = l1, l2 = l2, m1 = m1, m2 = m2, g = g, _X_ = _X_, _U_ = _U_)
    return PB.SafetyProblem(sys, _I_, _S_, PB.Infinity())
end

function optimal_control_problem(; l1 = 1.0, l2 = 1.0, m1 = 1.0, m2 = 1.0, g = 9.81)
    _X_ = UT.HyperRectangle(
        SVector(-π/4.0, -π/4.0, -5.0, -5.0),
        SVector(π/4.0, π/4.0, 5.0, 5.0),
    )
    _U_ = UT.LazySetMinus(
        UT.HyperRectangle(SVector(-5.5), SVector(5.5)),
        UT.HyperRectangle(SVector(-1.5), SVector(1.5));
        _I_ = UT.HyperRectangle(
            SVector(-3.0*pi/180.0, -3.0*pi/180.0, -0.5, -0.5),
            SVector(3.0*pi/180.0, 3.0*pi/180.0, 0.5, 0.5),
        ),
        _T_ = UT.HyperRectangle(
            SVector(35.0*pi/180.0, -15.0*pi/180.0, -1.0, -1.0),
            SVector(45.0*pi/180.0, 15.0*pi/180.0, 1.0, 1.0),
        ),
    )
    _X_ = UT.HyperRectangle(SVector(-π/2.0, -π, -5.0, -5.0), SVector(π/2.0, π, 5.0, 5.0))
    _U_ = UT.LazySetMinus(
        UT.HyperRectangle(SVector(-5.5), SVector(5.5)),
        UT.HyperRectangle(SVector(-0.5), SVector(0.5)),
    )
    _I_ = UT.HyperRectangle(
        SVector(-3.0*pi/180.0, -3.0*pi/180.0, -0.5, -0.5),
        SVector(3.0*pi/180.0, 3.0*pi/180.0, 0.5, 0.5),
    )
    _T_ = UT.HyperRectangle(
        SVector(-π/2.0, pi-50.0*pi/180.0, -4.5, -4.5),
        SVector(π/2.0, pi+50.0*pi/180.0, 4.5, 4.5),
    )
    sys = system(; l1 = l1, l2 = l2, m1 = m1, m2 = m2, g = g, _X_ = _X_, _U_ = _U_)
    return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
end

end
