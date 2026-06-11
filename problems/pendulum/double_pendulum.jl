module DoublePendulum

using StaticArrays
using MathematicalSystems
using Dionysos
import Symbolics
import IntervalArithmetic as IA
const UT = Dionysos.Utils
const ST = Dionysos.System
const PB = Dionysos.Problem

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

function jacobian_bound(; l1 = 1.0, l2 = 1.0, m1 = 1.0, m2 = 1.0, g = 9.81, ωmax = 5.0) # c'est du chat gpt mais de tout façon c'est utilisé nul part
    αmin = max(m1, 1.0e-6)
    M = m1 + m2
    gloc = abs(g)

    cθ1 =
        (2.0 * M * gloc) / (max(l1, 1.0e-6) * αmin) + (2.0 * m2 * (l1 + l2) * ωmax^2) / αmin
    cθ2 =
        (2.0 * M * gloc) / (max(l2, 1.0e-6) * αmin) + (2.0 * m2 * (l1 + l2) * ωmax^2) / αmin
    cω1 = 1.0 + (2.0 * m2 * l1 * ωmax) / αmin
    cω2 = 1.0 + (2.0 * m2 * l2 * ωmax) / αmin

    return u -> @SMatrix [
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
        cθ1 cθ2 cω1 cω2
        cθ1 cθ2 cω1 cω2
    ]
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

function optimal_control_problem(;
    l1 = 1.0,
    l2 = 1.0,
    m1 = 1.0,
    m2 = 1.0,
    g = 9.81,
    objective = "reachability_up_nonconvex",
)
    if objective == "benchmark_up_convex" # j'ai ajouté un set, mais le problème me semble plus mainigfull
        _X_ = UT.HyperRectangle(SVector(-π, -π, -6.0, -6.0), SVector(π, π, 6.0, 6.0))
        _U_ = UT.HyperRectangle(SVector(-8.5), SVector(8.5))
        _I_ = UT.HyperRectangle(
            SVector(-10.0 * pi / 180.0, -10.0 * pi / 180.0, -0.5, -0.5),
            SVector(10.0 * pi / 180.0, 10.0 * pi / 180.0, 0.5, 0.5),
        )
        _T_ = UT.HyperRectangle(
            SVector(π - 25.0*pi/180.0, π - 25.0*pi/180.0, -5.0, -5.0),
            SVector(π + 25.0*pi/180.0, π + 25.0*pi/180.0, 5.0, 5.0),
        )
    else
        _X_ = UT.HyperRectangle(
            SVector(-π / 2.0, -π, -5.0, -5.0),
            SVector(π / 2.0, π, 5.0, 5.0),
        )
        _U_ = UT.LazySetMinus(
            UT.HyperRectangle(SVector(-5.5), SVector(5.5)),
            UT.HyperRectangle(SVector(-0.5), SVector(0.5)),
        )
        _I_ = UT.HyperRectangle(
            SVector(-3.0 * pi / 180.0, -3.0 * pi / 180.0, -0.5, -0.5),
            SVector(3.0 * pi / 180.0, 3.0 * pi / 180.0, 0.5, 0.5),
        )
        _T_ = UT.HyperRectangle(
            SVector(-π / 2.0, pi - 50.0 * pi / 180.0, -4.5, -4.5),
            SVector(π / 2.0, pi + 50.0 * pi / 180.0, 4.5, 4.5),
        )
    end

    sys = system(; l1 = l1, l2 = l2, m1 = m1, m2 = m2, g = g, _X_ = _X_, _U_ = _U_)
    return PB.OptimalControlProblem(sys, _I_, _T_, nothing, nothing, PB.Infinity())
end

function symbolic_system(
    _X_;
    l1 = 1.0,
    l2 = 1.0,
    m1 = 1.0,
    m2 = 1.0,
    g = 9.81,
    _U_ = UT.HyperRectangle(SVector(-5.5), SVector(5.5)),
    Ts::Float64 = 0.1,
    ΔX = IA.IntervalBox(IA.interval(-0.15, 0.15), 4),
    ΔU = IA.IntervalBox(IA.interval(-0.15, 0.15), 1),
    ΔW = IA.IntervalBox(IA.interval(0.0, 0.0), 1),
    rk4_num_substeps::Int = 1,
    obstacles = Any[],
)
    Symbolics.@variables θ1 θ2 ω1 ω2 τ w1 T
    x = [θ1; θ2; ω1; ω2]
    u = [τ]
    w = [w1]

    M = m1 + m2
    f_cont_expr(xloc, uloc) = begin
        Δθ = xloc[1] - xloc[2]
        α = m1 + m2 * sin(Δθ)^2
        [
            xloc[3]
            xloc[4]
            -sin(Δθ) * (m2 * l1 * xloc[3]^2 * cos(Δθ) + m2 * l2 * xloc[4]^2) -
            g * (M * sin(xloc[1]) - m2 * sin(xloc[2]) * cos(Δθ)) / (l1 * α) + uloc[1]
            sin(Δθ) * (M * l1 * xloc[3]^2 + m2 * l2 * xloc[4]^2 * cos(Δθ)) +
            g * (M * sin(xloc[1]) * cos(Δθ) - M * sin(xloc[2])) / (l2 * α)
        ]
    end

    f_disc = ST.runge_kutta4(f_cont_expr, x, u, T, rk4_num_substeps)

    fsymbolicT = eval(ST.build_function(f_disc, x, u, w, T)[1])
    fsymbolic = Symbolics.substitute(f_disc, Dict(T => Ts))

    Wset = UT.HyperRectangle(SVector(0.0), SVector(0.0))
    Uformat = UT.format_input_set(_U_)
    Wformat = UT.format_noise_set(Wset)

    f_cont_fun = dynamic(; l1 = l1, l2 = l2, m1 = m1, m2 = m2, g = g)
    function f_eval(xv, uv, _wv)
        xsv = SVector{4, Float64}(xv)
        usv = SVector{1, Float64}(uv)
        xnext = ST.runge_kutta4(f_cont_fun, xsv, usv, Ts, rk4_num_substeps)
        return collect(xnext)
    end
    function f_backward_eval(xv, uv, _wv)
        xsv = SVector{4, Float64}(xv)
        usv = SVector{1, Float64}(uv)
        xprev = ST.runge_kutta4(f_cont_fun, xsv, usv, -Ts, rk4_num_substeps)
        return collect(xprev)
    end

    return ST.SymbolicSystem(
        fsymbolicT,
        fsymbolic,
        Ts,
        length(x),
        length(u),
        length(w),
        x,
        u,
        w,
        ΔX,
        ΔU,
        ΔW,
        _X_,
        _U_,
        Wset,
        obstacles,
        f_eval,
        f_backward_eval,
        Uformat,
        Wformat,
    )
end

end
