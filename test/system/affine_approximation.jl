module TestMain
using Test

using Symbolics
import IntervalArithmetic as IA
using MathematicalSystems

using Dionysos
const DI = Dionysos
const ST = DI.System

const TOL = 1e-6

@testset "buildAffineApproximation" begin
    Symbolics.@variables xi, ui, wi
    x = [xi]
    u = [ui]
    w = [wi]
    f = [x[1] * x[1] + u[1] + w[1]]
    x̄ = [0.0]
    ū = [0.0]
    w̄ = [0.0]

    v = 10.0
    X = [(-v,), (v,)]
    U = [(-v,), (v,)]
    W = [(-v,), (v,)]

    approx_sys, L = ST.buildAffineApproximation(f, x, u, w, x̄, ū, w̄, X, U, W)

    @test isapprox(L, [2.0, 0.0, 0.0], atol = TOL)
    @test approx_sys.X === X
    @test approx_sys.U === U
    @test approx_sys.W === W
    @test isapprox(
        MathematicalSystems.successor(
            approx_sys,
            [0.0],
            [0.0],
            [0.0];
            check_constraints = false,
        ),
        [0.0],
        atol = TOL,
    )
end

end
