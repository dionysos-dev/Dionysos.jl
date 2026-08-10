# The certification chain: iterate steps backward from the terminal ellipsoid,
# apply the soundness gates, and assemble the result.

"""
    CertificationResult

Outcome of a certification chain.

- `success` — every step certified *and* every enabled gate passed (including
  terminal containment when checked);
- `failed_k` — first failing step index, or `nothing`;
- `steps` — forward-ordered [`StepRecord`](@ref)s;
- `controller` — the per-step affine controllers `κ_1..κ_K` on success, else `nothing`;
- `terminal_contained` — whether the terminal ellipsoid lies inside the problem's
  target set (`nothing` when the check is disabled);
- `initial_coverage` — max of `(v − c₁)ᵀP₁(v − c₁)` over the initial set's vertices
  (≤ 1 means the entry funnel covers the initial set; `nothing` if unavailable);
- `state_domain_checked` — whether the reach-avoid gate could run on this domain type;
- `lmi_data` — `(; ellipsoids, kappas)` with the forward-ordered funnel.
"""
struct CertificationResult{S, CTRL, LMI}
    success::Bool
    failed_k::Union{Nothing, Int}
    solve_time_sec::Float64
    steps::Vector{S}
    controller::CTRL
    terminal_contained::Union{Nothing, Bool}
    initial_coverage::Union{Nothing, Float64}
    state_domain_checked::Bool
    lmi_data::LMI
end

function _collect_kappas(steps::AbstractVector{<:StepRecord})
    valid_steps = filter(step -> step.kappa !== nothing, steps)
    return [step.kappa for step in valid_steps]
end

function _terminal_ellipsoid(ctx::ChainContext)
    opts = ctx.options
    nx = length(ctx.xs[end])
    xN = ctx.xs[end]

    if opts.terminal_shape !== nothing
        # `terminal_shape` is the LazySets shape matrix Q of the terminal
        # ellipsoid {x : (x−c)ᵀQ⁻¹(x−c) ≤ 1} (semi-axes² on the diagonal).
        terminal_shape = Matrix{Float64}(opts.terminal_shape)
        @assert size(terminal_shape) == (nx, nx) "terminal_shape must have size ($nx, $nx)."
        return LazySets.Ellipsoid(collect(float.(xN)), terminal_shape), nothing
    end

    return _default_terminal_ellipsoid(ctx.problem.target_set, xN, opts.terminal_shrink)
end

function _assemble_result(
    success,
    failed_k,
    t0,
    steps,
    ellipsoids,
    terminal_contained,
    initial_cov,
    domain_checked,
)
    steps_forward = reverse(steps)
    ellipsoids_forward = reverse(ellipsoids)
    kappas = _collect_kappas(steps_forward)

    return CertificationResult(
        success,
        failed_k,
        Float64(time() - t0),
        steps_forward,
        success ? kappas : nothing,
        terminal_contained,
        initial_cov,
        domain_checked,
        (; ellipsoids = ellipsoids_forward, kappas),
    )
end

function run_chain!(ctx::ChainContext)
    t0 = time()
    opts = ctx.options

    E_next, terminal_reason = _terminal_ellipsoid(ctx)
    if E_next === nothing
        return CertificationResult(
            false,
            ctx.K,
            Float64(time() - t0),
            StepRecord[],
            nothing,
            false,
            nothing,
            false,
            (; ellipsoids = [], kappas = [], terminal_reason),
        )
    end

    terminal_contained =
        opts.check_terminal ? terminal_containment(E_next, ctx.problem.target_set) : nothing

    domain_checked =
        opts.check_state_domain && _state_domain_supported(ctx.problem.system.X)

    steps = StepRecord[]
    ellipsoids = [E_next]

    for k in ctx.K:-1:1
        rec = backward_step!(ctx, k, E_next)
        rec = apply_gates(rec, ctx)
        push!(steps, rec)

        if rec.status != :ok
            return _assemble_result(
                false,
                k,
                t0,
                steps,
                ellipsoids,
                terminal_contained,
                nothing,
                domain_checked,
            )
        end

        E_next = rec.ellipsoid
        push!(ellipsoids, rec.ellipsoid)
    end

    initial_cov = initial_coverage(ctx.problem.initial_set, ellipsoids[end])

    success = terminal_contained !== false

    return _assemble_result(
        success,
        success ? nothing : ctx.K + 1,
        t0,
        steps,
        ellipsoids,
        terminal_contained,
        initial_cov,
        domain_checked,
    )
end
