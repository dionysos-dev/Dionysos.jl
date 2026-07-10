module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Dionysos
using Clarabel, JuMP
import MathematicalSystems as MS
import LazySets
using HybridSystems

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System

println("Started test")

opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

# Shared fixture. The expected values below are golden outputs captured from the
# pre-redesign kernels (transition_fixed / _has_transition / transition_backward)
# on this exact problem — the redesigned API must reproduce them.
A = [1.0 0.1; -0.05 0.98]
B = reshape([0.0; 0.1], 2, 1)
g = [0.01; -0.02]
D = [0.005 0.0; 0.0 0.005]

Uset = UT.box(SVector(-5.0), SVector(5.0))
Uformat = ST.format_input_set(Uset)
Wrect = UT.box(SVector(-1.0, -1.0), SVector(1.0, 1.0))
Wmat = ST.format_noise_set(Wrect)

cost_fun = UT.QuadraticStateControlFunction(
    Matrix(1.0 * I, 2, 2),
    Matrix(1.0 * I, 1, 1),
    zeros(2, 1),
    zeros(2, 1),
    zeros(1, 1),
    1.0,
)
Λ = UT.get_full_psd_matrix(cost_fun)

# Golden fixtures were defined in the quadratic-form convention (P); LazySets
# stores Q = P⁻¹, so the source shape is inverted here.
E1 = LazySets.Ellipsoid([1.0, 1.0], inv([2.0 0.2; 0.2 1.5]))
E2 = LazySets.Ellipsoid([1.08, 0.9], [1.0 0.0; 0.0 1.0])

sys_noisy =
    MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, g, D, nothing, nothing, nothing)
sys_plain = HybridSystems.ConstrainedAffineControlDiscreteSystem(A, B, g, nothing, nothing)

@testset "solve_transition (noisy system: W mapped through D)" begin
    result = ST.solve_transition(sys_noisy, E1, E2, Uformat, Wmat, Λ, opt_sdp)
    @test result isa ST.TransitionResult
    @test result.feasible
    @test result.source === nothing
    # The controller is one of several near-optimizers (the cost is pinned much
    # tighter): the vertex ordering of the noise polytope shifts which one the
    # SDP returns, so the controller goldens carry a looser tolerance.
    @test result.cost ≈ 5.593497169 rtol = 1e-5
    @test vec(Matrix(result.controller.A)) ≈ [0.1035123589, -0.6746765442] atol = 1e-2
    @test collect(result.controller.c) ≈ [0.9701680955] atol = 1e-2

    # cost passed as a QuadraticStateControlFunction is equivalent
    result_f = ST.solve_transition(sys_noisy, E1, E2, Uformat, Wmat, cost_fun, opt_sdp)
    @test result_f.cost ≈ result.cost rtol = 1e-8

    # input constraints passed as a set are equivalent
    result_set = ST.solve_transition(sys_noisy, E1, E2, Uset, Wrect, Λ, opt_sdp)
    @test result_set.cost ≈ result.cost rtol = 1e-8
end

@testset "solve_transition (plain system: W used as given)" begin
    Wsmall = 0.005 .* Wmat
    result = ST.solve_transition(sys_plain, E1, E2, Uformat, Wsmall, Λ, opt_sdp)
    @test result.feasible
    @test result.cost ≈ 5.59349717 rtol = 1e-5
    @test vec(Matrix(result.controller.A)) ≈ [0.1031962696, -0.674806461] atol = 1e-2
end

@testset "solve_transition (infeasible)" begin
    # target far outside the one-step reachable set
    E_far = LazySets.Ellipsoid([1e6, 1e6], [0.01 0.0; 0.0 0.01])
    result = ST.solve_transition(sys_noisy, E1, E_far, Uformat, Wmat, Λ, opt_sdp)
    @test !result.feasible
    @test result.controller === nothing
    @test result.cost === nothing
end

@testset "solve_transition_backward" begin
    Lip = [1.0, 1.0, 0.5, 0.1]
    result = ST.solve_transition_backward(
        sys_noisy,
        E2,
        [1.0, 1.0],
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        maxδx = 2.0,
        maxδu = 2.0,
        λ = 0.01,
    )
    @test result.feasible
    @test result.source isa LazySets.Ellipsoid
    @test collect(LazySets.center(result.source)) ≈ [1.0, 1.0]
    @test Matrix(UT.get_quadratic_form(result.source)) ≈ [3.41386039 0.0; 0.0 3.41386039] atol =
        1e-3
    @test result.cost ≈ 4.824287182 rtol = 1e-4

    # The gains are only loosely pinned by the SDP (flat objective direction),
    # so check the certified property instead of golden gain values: the
    # controller maps every sampled source state into the target, for every
    # noise vertex. The Lipschitz term makes the certificate conservative, so
    # no slack is needed.
    K = Matrix(result.controller.A)
    b = collect(result.controller.c)
    for x in UT.samples(result.source, 100)
        u = K * x + b
        for j in 1:size(Wmat, 2)
            xnext = A * x + B * u + g + D * Wmat[:, j]
            @test xnext ∈ E2
        end
    end
end

@testset "stabilizing_feedback" begin
    A3 = [
        0.0 1.0 0.0
        0.0 0.0 1.0
        2.0 1.0 5.0
    ]
    B3 = [0.0; 0.0; 1.0][:, :]
    c3 = zeros(3)

    sys = HybridSystems.ConstrainedAffineControlDiscreteSystem(A3, B3, c3, Nothing, Nothing)
    feasible, K, P, gamma = ST.stabilizing_feedback(sys, opt_sdp)
    @test feasible
    @test K ≈ [-2.0 -1.0 -5.0] atol = 1e-1
    @test isposdef(Symmetric(Matrix(P)))
end

println("End test")

end  # module TestMain
