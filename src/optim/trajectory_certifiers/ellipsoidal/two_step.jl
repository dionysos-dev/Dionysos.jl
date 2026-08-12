# Two-step rescue (the lever against the per-step wall on underactuated plants):
# when the one-step transition at k dies at the maximal box, retry as a SINGLE
# two-step transition into E_{k+2} through the already-synthesized κ_{k+1} —
# no containment requirement at the (possibly needle-thin) intermediate funnel.
# See `ST.solve_transition_backward_2step`. Step 2's linearization box is fixed
# in advance (ΔX_max/ΔU_max of the adaptive options); its consistency — the
# realized intermediate excursion staying inside it — is checked a posteriori,
# in the house adaptive-box contract.

# Per-axis excursion radius of the intermediate state x₁ = Ã₁x₀ + b + e₁ over
# the synthesized source (shape Q₁) plus the step-1 remainder box.
function _intermediate_excursion(Ã₁, Q1, Lip1, δ_scale, nx)
    S = Ã₁ * Q1 * transpose(Ã₁)
    return [sqrt(max(0.0, S[i, i])) + Lip1[i] * δ_scale for i in 1:nx]
end

function _two_step_rescue!(ctx::ChainContext, k::Int, steps, ellipsoids)
    opts = ctx.options
    aopts = opts.adaptive_boxes
    aopts === nothing && return nothing
    # The 2-step kernel is :vertices-only (documented on ChainOptions and on
    # `ST.solve_transition_backward_2step`); other remainder chains skip the
    # rescue instead of crashing mid-chain or being silently downgraded.
    opts.remainder_model === :vertices || return nothing
    k <= ctx.K - 1 || return nothing                   # needs κ_{k+1}
    length(ellipsoids) >= 2 || return nothing          # needs E_{k+2}

    κ1 = steps[end].kappa
    κ1 isa MS.AffineMap || return nothing
    E_target = ellipsoids[end - 1]
    E_mid = ellipsoids[end]

    xk = collect(Float64, ctx.xs[k])
    uk = collect(Float64, ctx.us[k])
    x1 = collect(Float64, ctx.xs[k + 1])
    u1 = collect(Float64, ctx.us[k + 1])
    nx = length(xk)

    K1 = Matrix{Float64}(κ1.A)
    cap0 = opts.domain_cap ? _domain_cap(ctx.problem.system.X, xk) : nothing

    # Step 2's box must be sized to the ACTUAL intermediate excursion, not the
    # adaptive maxima: e₂'s CONSTANT budget scales with (box₂ radius)² and its
    # Lipschitz bound grows with the box, so an oversized box₂ makes the kernel
    # trivially infeasible (measured: adaptive-maxima sizing gave e₂ half-widths
    # six orders of magnitude above the target radii). Ladder it in multiples of
    # the skipped funnel E_{k+1} — the excursion the rescue exists to allow.
    r_mid = sqrt(
        max(
            1e-12,
            maximum(
                LA.eigvals(LA.Symmetric(Matrix{Float64}(LazySets.shape_matrix(E_mid)))),
            ),
        ),
    )

    for s2 in (8.0, 32.0, 128.0), scale in (1.0, 4.0)
        δx2 = clamp.(fill(s2 * r_mid, nx), aopts.ΔX_min, aopts.ΔX_max)
        δu2 = clamp.(abs.(K1) * δx2, aopts.ΔU_min, aopts.ΔU_max)
        approx2 =
            ST.build_affine_approximation(ctx.affine_provider, x1, u1; δx = δx2, δu = δu2)
        Lip2 = collect(Float64, approx2.lipschitz[1:nx])
        e2_scale = maximum(δx2)^2 + maximum(δu2)^2
        e2_hw = Lip2 .* e2_scale
        δx1 = _clamp_box_radii(scale .* aopts.ΔX_initial, aopts.ΔX_min, aopts.ΔX_max)
        δu1 = _clamp_box_radii(scale .* aopts.ΔU_initial, aopts.ΔU_min, aopts.ΔU_max)
        approx1 =
            ST.build_affine_approximation(ctx.affine_provider, xk, uk; δx = δx1, δu = δu1)

        source_cap = cap0
        box_shaved = δx1 .* (1 - 1e-6)
        source_cap = source_cap === nothing ? box_shaved : min.(source_cap, box_shaved)

        result = ST.solve_transition_backward_2step(
            approx1.system,
            approx2.system,
            κ1,
            E_target,
            xk,
            uk,
            approx1.Uformat,
            approx1.Wformat,
            approx2.Wformat,
            ctx.S,
            approx1.lipschitz,
            e2_hw,
            ctx.backend;
            maxδx = opts.maxδx,
            maxδu = opts.maxδu,
            λ = opts.λ,
            objective = opts.objective,
            source_cap = source_cap,
        )
        result.feasible || continue

        E_prev = result.source
        kappa = result.controller
        Q1 = Matrix{Float64}(LazySets.shape_matrix(E_prev))

        # Step-1 box consistency (the source and its controller image must live
        # inside the box the step-1 Hessian bound was taken on).
        required_X, required_U = _required_linearization_box_radii(E_prev, kappa, xk, uk)
        _box_contains_required_radii(δx1, δu1, required_X, required_U; atol = aopts.atol) ||
            continue

        # Step-2 box consistency: the realized intermediate excursion (source
        # through the closed step-1 loop, plus the step-1 remainder) must stay
        # inside step 2's fixed state box, and κ₁'s image over it inside its
        # input box.
        K0 = Matrix{Float64}(kappa.A)
        Ã₁ = approx1.system.A + approx1.system.B * K0
        Lip1v = collect(Float64, approx1.lipschitz[1:nx])
        δ_scale = maximum(δx1)^2 + maximum(δu1)^2
        r1 = _intermediate_excursion(Ã₁, Q1, Lip1v, δ_scale, nx)
        all(r1 .<= δx2 .+ aopts.atol) || continue
        u_dev = abs.(K1) * r1
        all(u_dev .<= δu2 .+ aopts.atol) || continue

        return StepRecord(
            k,
            :ok,
            Float64(result.cost),
            E_prev,
            kappa,
            _step_summary(
                approx1.lipschitz,
                δx1,
                δu1,
                required_X,
                required_U,
                1,
                :two_step_rescue;
                provider_summary = approx1.summary,
                _candidate_diagnostics_empty()...,
            ),
        )
    end

    return nothing
end
