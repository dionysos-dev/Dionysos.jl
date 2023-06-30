module NonLinear

using StaticArrays, LinearAlgebra, Symbolics, IntervalArithmetic
using MathematicalSystems, HybridSystems

import Dionysos
using Dionysos.Control
using Dionysos.Problem
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic
const PR = DI.Problem

function unstableSimple()
    Symbolics.@variables px py vx vy wx wy T

    f = [
        1.1 * px - 0.2 * py - 0.00005 * py^3 + T * vx
        1.1 * py + 0.2 * px + 0.00005 * px^3 + T * vy
    ]

    x = [px; py] # state
    u = [vx; vy] # control
    w = [wx; wy]
    return f, x, u, w, T
end

function system(X, U, W, obstacles, Ts)
    f, x, u, w, T = unstableSimple()
    fsymbolicT = eval(build_function(f, x, u, w, T)[1])
    #### PWA approximation description #####
    fsymbolic = Symbolics.substitute(f, Dict([T => Ts]))
    ΔX = IntervalBox(-1.0 .. 1.0, 2)
    ΔU = IntervalBox(-10 * 2 .. 10 * 2, 2)
    ΔW = IntervalBox(-0.0 .. 0.0, 1)

    ########## Inputs description ##########
    Uformat = SY.convert_input_set_into_format(U)
    function f_eval(x, u, w)
        return [
            1.1 * x[1] - 0.2 * x[2] - 0.00005 * x[2]^3 + Ts * u[1]
            1.1 * x[2] + 0.2 * x[1] + 0.00005 * x[1]^3 + Ts * u[2]
        ]
    end

    function f_backward_eval(x, u, w)
        return [
            1.1 * x[1] - 0.2 * x[2] - 0.00005 * x[2]^3 - Ts * u[1]
            1.1 * x[2] + 0.2 * x[1] + 0.00005 * x[1]^3 - Ts * u[2]
        ]
    end

    # f_eval(x, u, w) = fsymbolicT(x, u, w, Ts)
    # f_backward_eval(x, u, w) = fsymbolicT(x, u, w, -Ts)

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
        X,
        U,
        Uformat,
        W,
        obstacles,
        f_eval,
        f_backward_eval,
    )
end

function problem(;
    X = IntervalBox(-20.0 .. 20.0, 2),
    obstacles = [UT.Ellipsoid(Matrix{Float64}(I(2)) * 1 / 50, [0.0; 0.0])],
    U = UT.HyperRectangle(SVector(-10.0, -10.0), SVector(10.0, 10.0)),
    E0 = UT.Ellipsoid(Matrix{Float64}(I(2)) * 10.0, [-10.0; -10.0]),
    Ef = UT.Ellipsoid(Matrix{Float64}(I(2)) * 1.0, [10.0; 10.0]),
    state_cost = UT.ZeroFunction(),
    transition_cost = UT.QuadraticStateControlFunction(
        Matrix{Float64}(I(2)),
        Matrix{Float64}(I(2)),
        zeros(2, 2),
        zeros(2),
        zeros(2),
        1.0,
    ),
    W = 0.0 * [
        -1 -1 1 1
        -1 1 -1 1
    ],
    Ts = 1.0,
    N = Infinity(),
)
    sys = system(X, U, W, obstacles, Ts)

    problem = OptimalControlProblem(sys, E0, Ef, state_cost, transition_cost, N)
    return problem
end

end
