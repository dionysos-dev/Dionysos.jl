module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using Symbolics
import IntervalArithmetic as IA

const TOL = 1e-6

@testset "SymbolicAffineApproximationProvider" begin
    Symbolics.@variables xi, ui, wi
    x = [xi]
    u = [ui]
    w = [wi]
    f = [x[1] * x[1] + u[1] + w[1]]

    provider = ST.SymbolicAffineApproximationProvider(f, x, u, w, [10.0], nothing, nothing)

    approx = ST.build_affine_approximation(provider, [0.0], [0.0]; δx = [10.0], δu = [10.0])

    # Hessian of x² + u + w is 2 in x and 0 in u/w over the ±10 boxes.
    @test isapprox(approx.lipschitz, [2.0, 0.0, 0.0], atol = TOL)
    @test approx.summary.δx == [10.0]
    @test isapprox(
        MathematicalSystems.successor(
            approx.system,
            [0.0],
            [0.0],
            [0.0];
            check_constraints = false,
        ),
        [0.0],
        atol = TOL,
    )

    # Linearization at xbar = 1: A = [2], c = f(1,0,0) − A·1 = 1 − 2 = −1.
    approx1 = ST.build_affine_approximation(provider, [1.0], [0.0]; δx = [1.0], δu = [1.0])
    @test isapprox(approx1.system.A, [2.0;;], atol = TOL)
    @test isapprox(approx1.system.c, [-1.0], atol = TOL)
end

@testset "AnalyticAffineApproximationProvider" begin
    # Same dynamics with hand-written derivatives: f(x, u, w) = x² + u + w.
    provider = ST.AnalyticAffineApproximationProvider(;
        A = (xb, ub, wb) -> [2.0 * xb[1];;],
        B = (xb, ub, wb) -> [1.0;;],
        E = (xb, ub, wb) -> [1.0;;],
        f = (xb, ub, wb) -> [xb[1]^2 + ub[1] + wb[1]],
        lipschitz = (xb, ub, wb, δx, δu) -> [2.0, 0.0, 0.0],
        nw = 1,
        ΔW = [10.0],
        Uformat = nothing,
        Wformat = nothing,
    )

    approx = ST.build_affine_approximation(provider, [1.0], [0.0]; δx = [1.0], δu = [1.0])

    @test isapprox(approx.system.A, [2.0;;], atol = TOL)
    @test isapprox(approx.system.c, [-1.0], atol = TOL)
    @test approx.lipschitz == [2.0, 0.0, 0.0]
    @test isapprox(
        MathematicalSystems.successor(
            approx.system,
            [1.0],
            [0.0],
            [0.0];
            check_constraints = false,
        ),
        [1.0],
        atol = TOL,
    )

    # Constant-vector Lipschitz bound and default zero noise matrix.
    provider2 = ST.AnalyticAffineApproximationProvider(;
        A = (xb, ub, wb) -> [2.0 * xb[1];;],
        B = (xb, ub, wb) -> [1.0;;],
        f = (xb, ub, wb) -> [xb[1]^2 + ub[1]],
        lipschitz = [2.0, 0.0, 0.0],
        nw = 1,
        Uformat = nothing,
        Wformat = nothing,
    )
    approx2 = ST.build_affine_approximation(provider2, [1.0], [0.0]; δx = [1.0], δu = [1.0])
    @test approx2.lipschitz == [2.0, 0.0, 0.0]
    @test all(iszero, approx2.system.D)
end

@testset "compute_jacobian_bound (automatic growth bound)" begin
    # ẋ₁ = x₂, ẋ₂ = -sin(x₁) + u₁ over x₁ ∈ [-π/3, π/3]:
    # |∂f₂/∂x₁| = |-cos(x₁)| ∈ [1/2, 1] → bound 1; diagonal entries are 0.
    F_sys(x, u) = [x[2], -sin(x[1]) + u[1]]
    X = UT.box([-π / 3, -1.0], [π / 3, 1.0])
    U = UT.box([-1.0], [1.0])
    sys = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(F_sys, 2, 1, X, U)

    jb = ST.compute_jacobian_bound(sys)
    M = jb([0.5])
    @test isapprox(M[1, 1], 0.0; atol = TOL)
    @test isapprox(M[1, 2], 1.0; atol = TOL)
    @test isapprox(M[2, 1], 1.0; atol = TOL)
    @test isapprox(M[2, 2], 0.0; atol = TOL)

    # Signed diagonal supremum: ẋ = -x on [1, 2] has ∂f/∂x ≡ -1, not |−1|.
    F_dec(x, u) = [-x[1] + 0.0 * u[1]]
    sysd = MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        F_dec,
        1,
        1,
        UT.box([1.0], [2.0]),
        U,
    )
    Md = ST.compute_jacobian_bound(sysd)([0.0])
    @test isapprox(Md[1, 1], -1.0; atol = TOL)

    # The automatic bound plugs straight into the growth-bound constructor.
    approx = ST.ContinuousTimeGrowthBound(sys)
    @test approx isa ST.ContinuousTimeGrowthBound
    out =
        ST.get_over_approximation_map(approx)(UT.box([-0.1, -0.1], [0.1, 0.1]), [0.5], 0.1)
    @test out isa UT.Box || out isa Dionysos.LazySets.Hyperrectangle
end

end
