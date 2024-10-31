module SimpleProblem

using StaticArrays, LinearAlgebra, Symbolics, IntervalArithmetic
using MathematicalSystems, HybridSystems

import Dionysos
using Dionysos.Problem
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem

struct SimpleSystem{N, T} <: ST.ControlSystem{N, T}
    X::Any
    U::Any
    tstep::Float64
    measnoise::SVector{N, T}
    periodic::Any
    periods::Any
    T0::Any
    sys_map::Any
    f_eval::Any
    f_backward::Any
end

function system(
    rectX,
    obstacles,
    rectU,
    Uobstacles,
    tstep,
    measnoise,
    periodic,
    periods,
    T0,
)
    function sys_map(x::SVector{N, T}, u, tstep) where {N, T}
        return x + tstep * u
    end
    function f_eval(x::SVector{N, T}, u) where {N, T}
        return sys_map(x, u, tstep)
    end
    function f_backward(x::SVector{N, T}, u) where {N, T}
        return sys_map(x, u, -tstep)
    end
    obs = UT.LazyUnionSetArray(obstacles)
    X = UT.LazySetMinus(rectX, obs)
    Uobs = UT.LazyUnionSetArray(Uobstacles)
    U = UT.LazySetMinus(rectU, Uobs)
    return SimpleSystem(
        X,
        U,
        tstep,
        measnoise,
        periodic,
        periods,
        T0,
        sys_map,
        f_eval,
        f_backward,
    )
end

function problem(;
    rectX = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0)),
    obstacles = [UT.HyperRectangle(SVector(15.0, 15.0), SVector(20.0, 20.0))],
    periodic = Int[],
    periods = [30.0, 30.0],
    T0 = [0.0, 0.0],
    rectU = UT.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0)),
    Uobstacles = [UT.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))],
    _I_ = UT.HyperRectangle(SVector(5.0, 5.0), SVector(6.0, 6.0)),
    _T_ = UT.HyperRectangle(SVector(25.0, 25.0), SVector(28.0, 28.0)),
    state_cost = UT.ZeroFunction(),
    transition_cost = UT.ConstantControlFunction(0.5),
    tstep = 0.8,
    measnoise = SVector(0.0, 0.0),
    N = Infinity(),
)
    sys =
        system(rectX, obstacles, rectU, Uobstacles, tstep, measnoise, periodic, periods, T0)
    problem = OptimalControlProblem(sys, _I_, _T_, state_cost, transition_cost, N)
    return problem
end

end
