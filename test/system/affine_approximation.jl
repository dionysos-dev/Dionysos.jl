module TestMain
using Test

using Symbolics
import IntervalArithmetic as IA
using MathematicalSystems

using Dionysos
const DI = Dionysos
const ST = DI.System
const UT = DI.Utils

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
