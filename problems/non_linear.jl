module NonLinear

using StaticArrays
import LinearAlgebra as LA

import Symbolics
import IntervalArithmetic as IA

import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

function unstableSimple(; μ = 0.00005, noise = false)
    Symbolics.@variables px py vx vy wx wy T
    if noise
        f = [
            1.1 * px - 0.2 * py - μ * py^3 + T * vx + wx
            1.1 * py + 0.2 * px + μ * px^3 + T * vy + wy
        ]
    else
        f = [
            1.1 * px - 0.2 * py - μ * py^3 + T * vx
            1.1 * py + 0.2 * px + μ * px^3 + T * vy
        ]
    end

    x = [px; py]
    u = [vx; vy]
    w = [wx; wy]
    return f, x, u, w, T
end

function system(X, U, W, obstacles, Ts, noise, μ)
    f, x, u, w, T = unstableSimple(; noise = noise, μ = μ)

    fsymbolicT = eval(Symbolics.build_function(f, x, u, w, T)[1])

    #### Local approximation domains ####
    ΔX = [IA.interval(-1.0, 1.0), IA.interval(-1.0, 1.0)]
    ΔU = [IA.interval(-20.0, 20.0), IA.interval(-20.0, 20.0)]
    ΔW = [IA.interval(0.0, 0.0), IA.interval(0.0, 0.0)]

    fsymbolic = Symbolics.substitute(f, Dict(T => Ts))

    #### Format of input and noise set ####
    Uformat = UT.format_input_set(U)
    Wformat = UT.format_noise_set(W)

    #### Forward and backward dynamics ####
    function f_eval(x, u, w)
        return [
            1.1 * x[1] - 0.2 * x[2] - μ * x[2]^3 + Ts * u[1] + w[1]
            1.1 * x[2] + 0.2 * x[1] + μ * x[1]^3 + Ts * u[2] + w[2]
        ]
    end

    function f_backward_eval(x, u, w)
        return [
            1.1 * x[1] - 0.2 * x[2] - μ * x[2]^3 - Ts * u[1] - w[1]
            1.1 * x[2] + 0.2 * x[1] + μ * x[1]^3 - Ts * u[2] - w[2]
        ]
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
        X,
        U,
        W,
        obstacles,
        f_eval,
        f_backward_eval,
        Uformat,
        Wformat,
    )
end

function problem(;
    X = UT.HyperRectangle(SVector(-20.0, -20.0), SVector(20.0, 20.0)),

    obstacles = [UT.Ellipsoid(Matrix{Float64}(LA.I, 2, 2) * (1 / 30), [0.0; 0.0])],

    U = UT.HyperRectangle(SVector(-10.0, -10.0), SVector(10.0, 10.0)),

    E0 = UT.Ellipsoid(Matrix{Float64}(LA.I, 2, 2) * 10.0, [-10.0; -10.0]),

    Ef = UT.Ellipsoid(Matrix{Float64}(LA.I, 2, 2) * 1.0, [10.0; 10.0]),

    state_cost = UT.ZeroFunction(),

    transition_cost = UT.QuadraticStateControlFunction(
        Matrix{Float64}(LA.I, 2, 2),
        Matrix{Float64}(LA.I, 2, 2),
        zeros(2, 2),
        zeros(2),
        zeros(2),
        1.0,
    ),

    W = UT.HyperRectangle(SVector(0.0, 0.0), SVector(0.0, 0.0)),
    Ts = 1.0,
    noise = false,
    μ = 0.00005,
)
    sys = system(X, U, W, obstacles, Ts, noise, μ)
    return PR.OptimalControlProblem(sys, E0, Ef, state_cost, transition_cost)
end

end # module
