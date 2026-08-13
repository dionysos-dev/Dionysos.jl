# Two-step backward kernel: certify E_source → (κ₀, then the FIXED κ₁) → target
# in ONE LMI, with NO containment requirement at the intermediate state. This is
# the lever against the per-step bottleneck on rank-deficient-B plants: a
# per-step chain forces the one-step image inside the (possibly needle-thin)
# intermediate funnel, while the two-step transition only needs the TWO-step
# image inside the target — the certificate gains the 2-step coupling channel
# u → (coupled states). Convexity comes from κ₁ being data (the backward chain
# has already synthesized it): with M = A₂ + B₂K₁ constant, the composed system
# Ā = M·A₁, B̄ = M·B₁, c̄ = M·c₁ + B₂b₁ + c₂ stays linear in (L, F).
#
# Uncertainty bookkeeping: step-1's remainder e₁ (box Lip₁·(δx+δu), variable
# scale) passes through M, step-2's remainder e₂ is a CONSTANT box (its
# linearization box is fixed by the caller, who must verify a posteriori that
# the realized intermediate excursion stays inside it — the same consistency
# contract as the adaptive boxes). Aligned boxes sum: the combined vertex set is
# the usual 2ⁿ corners with per-axis half-widths |M|·Lip₁·(δx+δu) + e₂hw. The
# second input constraint |U·u₁| ≤ 1 is enforced over the true intermediate
# reach (u₁ = K₁x₁ + b₁ with x₁ affine in ξ and e₁): one block per (Uᵢ, e₁
# vertex). The cost block covers (x₀, u₀) only — the two-step transition cost
# of u₁ is not charged (documented limitation).

function _two_step_kernel(
    A1,
    B1,
    c1v,
    A2,
    B2,
    c2v,
    K1,
    b1,
    Wcols1,
    Wcols2,
    U,
    Λ,
    c0,
    u_ref0,
    cT,
    Q2,
    Lip1,
    e2_hw,
    sdp_solver;
    maxδx,
    maxδu,
    λ,
    objective::Symbol,
    source_cap = nothing,
)
    nx = length(c1v)
    nu = size(U[1], 2)
    size(Wcols1, 2) == 1 && size(Wcols2, 2) == 1 || error(
        "the two-step kernel currently supports a single noise vertex per step " *
        "(deterministic or fixed-offset noise).",
    )

    M = A2 + B2 * K1
    Ā = M * A1
    B̄ = M * B1
    c̄ = M * c1v + B2 * vec(b1) + c2v
    ν = vec(M * Wcols1[:, 1] + Wcols2[:, 1])

    # Combined remainder box: variable part |M|·Lip₁ (scaled by δx+δu), constant
    # part e₂hw. Corners are enumerated on the unit box once.
    hvar = abs.(M) * collect(Float64, Lip1[1:nx])
    hconst = collect(Float64, e2_hw)
    σs = collect.(LazySets.vertices_list(LazySets.Hyperrectangle(zeros(nx), ones(nx))))
    Nx = length(σs)
    Nu = length(U)

    km = _KernelModel(sdp_solver, _PSD_STRICTNESS)
    model = km.model
    @variable(model, L[i = 1:nx, j = 1:nx], PSD)
    @variable(model, F[i = 1:nu, j = 1:nx])
    @variable(model, ell[i = 1:nu, j = 1:1])
    @variable(model, beta[i = 1:Nx] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, tau2[i = 1:Nu, v = 1:Nx] >= 0)
    @variable(model, δx >= 0)
    @variable(model, δu >= 0)
    @variable(model, ϕ >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)
    _require_nonneg!(km, beta, tau, tau2)

    _pose_source_cap!(km, source_cap, L, nx)

    @expressions(model, begin
        At, Ā * L + B̄ * F
        ct, c̄ + B̄ * ell
    end)

    for (i, σ) in enumerate(σs)
        aux = @expression(
            model,
            Ā * hcat(c0) + hcat(ct) - hcat(cT) +
            hcat(σ .* hvar) * (δx + δu) +
            hcat(σ .* hconst .+ ν)
        )
        _pose_psd!(km, _reach_block, beta[i], _eye(nx), At, aux, Q2, nx)
    end

    # Input feasibility of u₀ over the source.
    _pose_inputs!(km, tau, U, _eye(nx), F, ell, nx)

    # Input feasibility of u₁ = K₁x₁ + b₁ over the intermediate reach:
    # x₁ = A₁x₀ + B₁u₀ + c₁ + e₁, so the ξ-gain is K₁(A₁L + B₁F) and the offset
    # carries the nominal-with-ℓ intermediate state plus each e₁ corner.
    K1AL = @expression(model, K1 * (A1 * L + B1 * F))
    x1ℓ = @expression(model, A1 * hcat(c0) + B1 * ell + hcat(c1v))
    Lip1v = collect(Float64, Lip1[1:nx])
    for i in 1:Nu, (v, σ) in enumerate(σs)
        Uoff = @expression(
            model,
            U[i] * (K1 * x1ℓ + hcat(vec(b1)) + hcat(K1 * (σ .* Lip1v)) * (δx + δu))
        )
        _pose_psd!(km, _input_block, tau2[i, v], _eye(nx), U[i] * K1AL, Uoff, nx)
    end

    z = zeros(nx, 1)
    G = [transpose(L) transpose(F) z]
    d = [transpose(c0) transpose(ell) 1]
    _pose_cost!(km, γ, J, _eye(nx), G, d, Λ, nx)

    _pose_proximity!(km, ϕ, _eye(nx), F, ell, u_ref0, δu, nx, nu)
    _pose_psd!(km, _source_radius_block, L, δx, nx)
    _pose_deviation_caps!(km, δx, maxδx, δu, maxδu)

    _pose_source_size_objective!(km, objective, λ, L, J, nx)

    optimize!(model)

    if !_solver_accepted(model)
        @debug "solve_transition_backward_2step: infeasible SDP" termination_status(model)
        return false, nothing, nothing, nothing
    end
    _validated(km, "solve_transition_backward_2step") ||
        return false, nothing, nothing, nothing

    Q1, kappa = _source_solution(L, F, ell, "solve_transition_backward_2step")
    Q1 === nothing && return false, nothing, nothing, nothing
    return true, Q1, kappa, value(J)
end

"""
    solve_transition_backward_2step(affsys1, affsys2, κ1, target, source_center,
                                    u_ref0, U, W1, W2, cost, lipschitz1,
                                    e2_halfwidths, sdp_solver;
                                    maxδx = 100.0, maxδu = 20.0, λ = 0.01,
                                    objective = :maximin,
                                    remainder_model = :vertices,
                                    source_cap = nothing)
        -> TransitionResult

Synthesize a first-step controller `κ₀` *and the largest source ellipsoid*
centered at `source_center` that reaches `target` in TWO steps: one step of
`affsys1` under `κ₀`, then one step of `affsys2` under the **fixed** controller
`κ1` (a `MathematicalSystems.AffineMap`, `u₁ = K₁x₁ + b₁`). No containment is
required at the intermediate state — the lever against per-step chain
bottlenecks on underactuated plants (the certificate gains the two-step
coupling channel).

`lipschitz1` bounds step 1's linearization error (scaled by the synthesized
`δx + δu` — SQUARED deviations, as in [`solve_transition_backward`](@ref));
`e2_halfwidths` is the CONSTANT per-axis half-width vector of step 2's
remainder box — the caller fixes step 2's linearization box in advance, scales
`e2_halfwidths` with the SQUARED box radii, and must verify a posteriori that
the realized intermediate excursion stays inside the box. Both steps currently
require a single noise vertex, and only `remainder_model = :vertices` is
supported (the combined two-box remainder has no ball form yet — the kwarg
exists so callers passing another model fail loudly instead of being silently
downgraded). The transition-cost bound covers `(x₀, u₀)` only. Returns a
[`TransitionResult`](@ref) whose `controller` is `κ₀` and whose `source` is the
synthesized two-step funnel entry.
"""
function solve_transition_backward_2step(
    affsys1::AffineSys,
    affsys2::AffineSys,
    κ1::MS.AffineMap,
    target::LazySets.Ellipsoid,
    source_center,
    u_ref0,
    U,
    W1,
    W2,
    cost,
    lipschitz1,
    e2_halfwidths,
    sdp_solver;
    maxδx = 100.0,
    maxδu = 20.0,
    λ = 0.01,
    objective::Symbol = :maximin,
    remainder_model::Symbol = :vertices,
    source_cap = nothing,
)
    remainder_model === :vertices || error(
        "solve_transition_backward_2step supports remainder_model = :vertices " *
        "only (the combined two-box remainder has no ball form), " *
        "got $remainder_model.",
    )
    feasible, Q1, kappa, J = _two_step_kernel(
        affsys1.A,
        affsys1.B,
        affsys1.c,
        affsys2.A,
        affsys2.B,
        affsys2.c,
        Matrix(κ1.A),
        collect(κ1.c),
        _check_noise_vertices(_state_noise(affsys1, W1)),
        _check_noise_vertices(_state_noise(affsys2, W2)),
        _input_matrices(U),
        _cost_factor(cost),
        source_center,
        u_ref0,
        LazySets.center(target),
        LazySets.shape_matrix(target),
        lipschitz1,
        e2_halfwidths,
        sdp_solver;
        maxδx = maxδx,
        maxδu = maxδu,
        λ = λ,
        objective = objective,
        source_cap = source_cap,
    )
    feasible || return _infeasible_transition()
    K, ℓ = _split_controller(kappa)
    controller = MS.AffineMap(K, ℓ - K * source_center)
    source = LazySets.Ellipsoid(
        collect(float.(source_center)),
        UT._symmetrize(Q1);
        check_posdef = false,
    )
    return TransitionResult(true, controller, J, source)
end
