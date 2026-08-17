module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

using LinearAlgebra
using Clarabel, JuMP
import MathematicalSystems as MS
import LazySets
using HybridSystems

opt_sdp = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

# Shared fixture. The expected values below are golden outputs captured from the
# pre-redesign kernels (transition_fixed / _has_transition / transition_backward)
# on this exact problem — the redesigned API must reproduce them.
A = [1.0 0.1; -0.05 0.98]
B = reshape([0.0; 0.1], 2, 1)
g = [0.01; -0.02]
D = [0.005 0.0; 0.0 0.005]

Uset = LazySets.Hyperrectangle(; low = SVector(-5.0), high = SVector(5.0))
Uformat = ST.format_input_set(Uset)
Wrect = LazySets.Hyperrectangle(; low = SVector(-1.0, -1.0), high = SVector(1.0, 1.0))
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
    @test result.cost ≈ 5.593497169 rtol = 1e-5
    # The controller is one of several near-optimizers (only the cost is pinned
    # tight), and which one the SDP returns shifts with the platform and the
    # PSD strictness margin — so check the certified property instead of golden
    # gain values: the controller maps every sampled source state into the
    # target, for every noise vertex.
    K = Matrix(result.controller.A)
    b = collect(result.controller.c)
    for x in UT.samples(E1, 100)
        u = K * x + b
        for j in 1:size(Wmat, 2)
            @test A * x + B * u + g + D * Wmat[:, j] ∈ E2
        end
    end

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
    # Certified property in place of golden gains (see the noisy testset).
    K = Matrix(result.controller.A)
    b = collect(result.controller.c)
    for x in UT.samples(E1, 100)
        u = K * x + b
        for j in 1:size(Wsmall, 2)
            @test A * x + B * u + g + Wsmall[:, j] ∈ E2
        end
    end
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

@testset "solve_transition_backward (:maximin objective)" begin
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
        objective = :maximin,
    )
    @test result.feasible
    Q = Matrix(LazySets.shape_matrix(result.source))
    evals = eigvals(Symmetric(Q))
    # The whole point of :maximin — no pancake: the smallest semi-axis is genuinely
    # positive and the aspect ratio stays modest.
    @test sqrt(minimum(evals)) > 1e-3
    @test maximum(evals) / minimum(evals) < 100.0

    # Same certified property as the :logdet golden testset.
    K = Matrix(result.controller.A)
    b = collect(result.controller.c)
    for x in UT.samples(result.source, 50)
        u = K * x + b
        for j in 1:size(Wmat, 2)
            @test A * x + B * u + g + D * Wmat[:, j] ∈ E2
        end
    end
end

@testset "solve_transition_backward (source_cap slab)" begin
    Lip = [1.0, 1.0, 0.5, 0.1]
    common = (; maxδx = 2.0, maxδu = 2.0, λ = 0.01, objective = :maximin)
    free = ST.solve_transition_backward(
        sys_noisy,
        E2,
        [1.0, 1.0],
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        common...,
    )
    @test free.feasible
    # Cap the second axis well below the uncapped radius so the constraint binds.
    r_free = sqrt.(diag(Matrix(LazySets.shape_matrix(free.source))))
    cap = [r_free[1] * 2.0, r_free[2] / 2.0]
    capped = ST.solve_transition_backward(
        sys_noisy,
        E2,
        [1.0, 1.0],
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        source_cap = cap,
        common...,
    )
    @test capped.feasible
    Qc = Matrix(LazySets.shape_matrix(capped.source))
    # Q[i,i] ≤ dᵢ² is exactly "ellipsoid ⊆ slab |xᵢ − c₁ᵢ| ≤ dᵢ".
    @test all(sqrt.(diag(Qc)) .<= cap .+ 1e-6)

    # The capped source still certifies the transition.
    K = Matrix(capped.controller.A)
    b = collect(capped.controller.c)
    for x in UT.samples(capped.source, 50)
        u = K * x + b
        for j in 1:size(Wmat, 2)
            @test A * x + B * u + g + D * Wmat[:, j] ∈ E2
        end
    end
end

@testset "solve_transition_forward (duality with backward)" begin
    Lip = [1.0, 1.0, 0.5, 0.1]
    result_b = ST.solve_transition_backward(
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
    @test result_b.feasible

    # Duality: backward certified that its synthesized source reaches E2, and its
    # remainder bound dominates forward's (δx_const = λ_max(Q₁) ≤ backward's δx
    # variable) — so the forward SDP from that source, with E2's shape, must be
    # feasible with contraction α ≤ 1 (mod solver tolerance). This cross-validates
    # both kernels, both remainder models, and the shared S-procedure blocks.
    Qhat = Matrix(LazySets.shape_matrix(E2))
    result_f = ST.solve_transition_forward(
        sys_noisy,
        result_b.source,
        collect(LazySets.center(E2)),
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        target_shape = Qhat,
        maxδu = 2.0,
    )
    @test result_f.feasible
    @test result_f.target isa LazySets.Ellipsoid
    α = tr(Matrix(LazySets.shape_matrix(result_f.target))) / tr(Qhat)
    @test α <= 1.0 + 1e-6

    # Certified property, forward direction: every sampled source state lands in
    # the synthesized target under the forward controller, for every noise vertex.
    K = Matrix(result_f.controller.A)
    b = collect(result_f.controller.c)
    for x in UT.samples(result_b.source, 50)
        u = K * x + b
        for j in 1:size(Wmat, 2)
            @test A * x + B * u + g + D * Wmat[:, j] ∈ result_f.target
        end
    end

    # Free-shape mode: tighter or equal (in trace) than the α-scaled fixed shape,
    # inside the conditioning sandwich.
    result_free = ST.solve_transition_forward(
        sys_noisy,
        result_b.source,
        collect(LazySets.center(E2)),
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        maxδu = 2.0,
        q_min = 1e-6,
        q_max = 1e6,
    )
    @test result_free.feasible
    @test tr(Matrix(LazySets.shape_matrix(result_free.target))) <=
          tr(Matrix(LazySets.shape_matrix(result_f.target))) + 1e-6
end

@testset "remainder_model = :ball (backward + forward)" begin
    Lip = [1.0, 1.0, 0.5, 0.1]
    rb = ST.solve_transition_backward(
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
        remainder_model = :ball,
    )
    @test rb.feasible
    # Same certified property as the :vertices golden testset: the ball model is a
    # sound over-approximation, so every sampled source state must land in the
    # target under every noise vertex.
    K = Matrix(rb.controller.A)
    b = collect(rb.controller.c)
    for x in UT.samples(rb.source, 50)
        u = K * x + b
        for j in 1:size(Wmat, 2)
            @test A * x + B * u + g + D * Wmat[:, j] ∈ E2
        end
    end

    # :john_ball — the box's John ellipsoid (per-axis radii √n·Lipᵢ) instead of
    # the scalar ball: a sound cover, so the same certified property must hold.
    rj = ST.solve_transition_backward(
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
        remainder_model = :john_ball,
    )
    @test rj.feasible
    Kj = Matrix(rj.controller.A)
    bj = collect(rj.controller.c)
    for x in UT.samples(rj.source, 50)
        u = Kj * x + bj
        for j in 1:size(Wmat, 2)
            @test A * x + B * u + g + D * Wmat[:, j] ∈ E2
        end
    end
    # The forward kernel does not implement :john_ball — it must say so.
    @test_throws ErrorException ST.solve_transition_forward(
        sys_noisy,
        rb.source,
        collect(LazySets.center(E2)),
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        target_shape = Matrix(LazySets.shape_matrix(E2)),
        maxδu = 2.0,
        remainder_model = :john_ball,
    )

    # Duality holds under the ball model too (forward's constant δx is dominated
    # by backward's δx variable, same remainder model on both sides).
    rf = ST.solve_transition_forward(
        sys_noisy,
        rb.source,
        collect(LazySets.center(E2)),
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        target_shape = Matrix(LazySets.shape_matrix(E2)),
        maxδu = 2.0,
        remainder_model = :ball,
    )
    @test rf.feasible
    α = tr(Matrix(LazySets.shape_matrix(rf.target))) / tr(Matrix(LazySets.shape_matrix(E2)))
    @test α <= 1.0 + 1e-6
end

@testset "solve_transition_backward_2step" begin
    Lip = [1.0, 1.0, 0.5, 0.1]
    # The fixed second-step controller comes from a certified one-step backward
    # transition into E2.
    Wz = zeros(2, 1)     # the two-step kernel requires a single noise vertex
    rb1 = ST.solve_transition_backward(
        sys_noisy,
        E2,
        [1.0, 1.0],
        [0.0],
        Uformat,
        Wz,
        Λ,
        Lip,
        opt_sdp;
        maxδx = 2.0,
        maxδu = 2.0,
        λ = 0.01,
        objective = :maximin,
    )
    @test rb1.feasible

    # Two-step transition from a DIFFERENT center, through κ₁, into E2 — with no
    # containment requirement at the intermediate state.
    e2_hw = 0.2 .* Lip[1:2]
    r2 = ST.solve_transition_backward_2step(
        sys_noisy,
        sys_noisy,
        rb1.controller,
        E2,
        [0.9, 1.1],
        [0.0],
        Uformat,
        Wz,
        Wz,
        Λ,
        Lip,
        e2_hw,
        opt_sdp;
        maxδx = 2.0,
        maxδu = 2.0,
        λ = 0.01,
        objective = :maximin,
    )
    @test r2.feasible
    @test r2.source isa LazySets.Ellipsoid
    @test collect(LazySets.center(r2.source)) ≈ [0.9, 1.1]

    # Certified property, sampled: apply κ₀ then the fixed κ₁ through the exact
    # affine plant (zero linearization error, zero noise — the certificate must
    # hold a fortiori) and land in E2, with both inputs feasible.
    K0 = Matrix(r2.controller.A)
    b0 = collect(r2.controller.c)
    K1 = Matrix(rb1.controller.A)
    b1 = collect(rb1.controller.c)
    for x0 in UT.samples(r2.source, 50)
        u0 = K0 * x0 + b0
        @test all(abs.(u0) .<= 5.0 + 1e-6)
        x1 = A * x0 + B * u0 + g
        u1 = K1 * x1 + b1
        @test all(abs.(u1) .<= 5.0 + 1e-6)
        x2 = A * x1 + B * u1 + g
        @test x2 ∈ E2
    end

    # Multi-vertex noise is rejected explicitly.
    @test_throws ErrorException ST.solve_transition_backward_2step(
        sys_noisy,
        sys_noisy,
        rb1.controller,
        E2,
        [0.9, 1.1],
        [0.0],
        Uformat,
        Wmat,
        Wz,
        Λ,
        Lip,
        e2_hw,
        opt_sdp,
    )
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

    # Unsuccessful solve (iteration-starved solver): NO values may be extracted
    # from the failed model — the old code inverted `value.(S)` before checking.
    opt_starved =
        optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true, "max_iter" => 1)
    feas_bad, K_bad, P_bad, γ_bad = ST.stabilizing_feedback(sys, opt_starved)
    @test !feas_bad
    @test K_bad === nothing && P_bad === nothing && γ_bad === nothing
end

@testset "argument validation" begin
    x0 = [1.0, 1.0]
    Lip = [0.1, 0.1, 0.0, 0.0]

    # Unknown remainder models must error loudly, never silently fall back.
    @test_throws ErrorException ST.solve_transition_backward(
        sys_noisy,
        E2,
        x0,
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        remainder_model = :box,
    )
    @test_throws ErrorException ST.solve_transition_forward(
        sys_noisy,
        E1,
        [1.08, 0.9],
        [0.0],
        Uformat,
        Wmat,
        Λ,
        Lip,
        opt_sdp;
        remainder_model = :Ball,
    )

    # An empty disturbance-vertex set would drop every reachability constraint.
    W_empty = zeros(2, 0)
    @test_throws ErrorException ST.solve_transition(
        sys_noisy,
        E1,
        E2,
        Uformat,
        W_empty,
        Λ,
        opt_sdp,
    )

    # The two-step kernel is :vertices-only and says so.
    κ1 = MS.AffineMap([0.0 0.0], [0.0])
    @test_throws ErrorException ST.solve_transition_backward_2step(
        sys_plain,
        sys_plain,
        κ1,
        E2,
        x0,
        [0.0],
        Uformat,
        zeros(2, 1),
        zeros(2, 1),
        Λ,
        Lip,
        [0.01, 0.01],
        opt_sdp;
        remainder_model = :ball,
    )
end

end  # module TestMain
