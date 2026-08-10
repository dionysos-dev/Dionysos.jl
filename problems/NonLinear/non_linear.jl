module NonLinear

using StaticArrays
import LinearAlgebra as LA
import LazySets

import Symbolics

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

function system(X, U, W, Ts, noise, μ)
    f, x, u, w, T = unstableSimple(; noise = noise, μ = μ)

    #### Local approximation radii ####
    ΔX = [1.0, 1.0]
    ΔU = [20.0, 20.0]
    ΔW = [0.0, 0.0]

    fsymbolic = Symbolics.substitute(f, Dict(T => Ts))

    #### Format of input and noise set ####
    Uformat = ST.format_input_set(U)
    Wformat = ST.format_noise_set(W)

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
        fsymbolic,
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
        f_eval,
        f_backward_eval,
        Uformat,
        Wformat,
    )
end

# Default avoid set for this benchmark. Obstacles are not part of the problem
# specification types; pass them to the solver (e.g. the lazy-ellipsoids
# `obstacles` attribute) and to the plots.
default_obstacles() = [LazySets.Ellipsoid([0.0; 0.0], Matrix{Float64}(LA.I, 2, 2) * 30.0)]

function problem(;
    X = LazySets.Hyperrectangle(; low = SVector(-20.0, -20.0), high = SVector(20.0, 20.0)),

    U = LazySets.Hyperrectangle(; low = SVector(-10.0, -10.0), high = SVector(10.0, 10.0)),

    E0 = LazySets.Ellipsoid([-10.0; -10.0], Matrix{Float64}(LA.I, 2, 2) * 0.1),

    Ef = LazySets.Ellipsoid([10.0; 10.0], Matrix{Float64}(LA.I, 2, 2) * 1.0),

    state_cost = UT.ZeroFunction(),

    transition_cost = UT.QuadraticStateControlFunction(
        Matrix{Float64}(LA.I, 2, 2),
        Matrix{Float64}(LA.I, 2, 2),
        zeros(2, 2),
        zeros(2),
        zeros(2),
        1.0,
    ),

    W = LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(0.0, 0.0)),
    Ts = 1.0,
    noise = false,
    μ = 0.00005,
)
    sys = system(X, U, W, Ts, noise, μ)
    return PR.OptimalControlProblem(sys, E0, Ef, state_cost, transition_cost)
end

end # module
