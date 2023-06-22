module TestMain
using Test

using Symbolics
using IntervalArithmetic
using MathematicalSystems

using Dionysos
const DI = Dionysos
const ST = DI.System

const TOL = 1e-6

@testset "get_f_eval" begin
    struct SystemExample{N, T, F <: Function} <: ST.ControlSystem{N, T}
        X::Any
        U::Any
        f_eval::F
    end
    function SystemExample(f::F) where {F <: Function}
        return SystemExample{1, Float64, F}((-Inf, Inf), (-Inf, Inf), f)
    end
    f(x, u) = 2 * x + u
    sys = SystemExample(f)
    x = 0.1
    u = 0.3
    @test (ST.get_f_eval)(sys)(x, u) === f(x, u)
end

@testset "BoundSecondOrder (a ≈ 0)" begin
    a = 0
    b = 1
    tstep = 0.1
    @test ST.BoundSecondOrder(a, b, tstep) === b * tstep
end

@testset "buildAffineApproximation" begin
    Symbolics.@variables xi, ui, wi
    x = [xi]
    u = [ui]
    w = [wi]
    f = [x[1] * x[1] + u[1] + w[1]]
    x̄ = [0.0]
    ū = [0.0]
    w̄ = [0.0]

    X = IntervalBox(-10.0 .. 10.0)
    U = IntervalBox(-10.0 .. 10.0)
    W = IntervalBox(-10.0 .. 10.0)
    approx_sys, L = ST.buildAffineApproximation(f, x, u, w, x̄, ū, w̄, X, U, W)

    @test isapprox(L, [2.0, 0.0, 0.0], atol = TOL)
    @test approx_sys.X === X
    @test approx_sys.U === U
    @test approx_sys.W === W
    @test isapprox(
        MathematicalSystems.successor(approx_sys, [0.0], [0.0], [0.0]),
        [0.0],
        atol = TOL,
    )
end

@testset "AffineApproximationDiscreteSystem" begin
    A = [2.0 3.0; 1.0 0.0]
    B = [1.0 0.0; 0.0 1.0]
    D = [1.0 0.0; 0.0 1.0]
    c = [1.0; 1.0]
    X = [(-Inf, -Inf), (Inf, Inf)]
    U = [(-Inf, -Inf), (Inf, Inf)]
    W = [(-Inf, -Inf), (Inf, Inf)]
    sys =
        MathematicalSystems.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, D, X, U, W)
    L = 1.0 # Just for the example, not mathematically correct

    obj_from_sys = ST.AffineApproximationDiscreteSystem(sys, L)
    @test obj_from_sys.constrainedAffineSys === sys
    @test obj_from_sys.L === L
    @test isapprox((obj_from_sys.f_eval)([0.0, 0.0], [0.0, 0.0], [0.0, 0.0]), c, atol = TOL)

    obj_from_matrices = ST.AffineApproximationDiscreteSystem(A, B, c, D, X, U, W, L)
    @test obj_from_matrices.constrainedAffineSys === sys
    @test obj_from_matrices.L === L
    @test isapprox(
        (obj_from_matrices.f_eval)([0.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
        c,
        atol = TOL,
    )
end

end
