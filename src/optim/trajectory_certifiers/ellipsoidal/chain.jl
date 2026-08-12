# The certification chain: iterate steps backward from the terminal ellipsoid,
# apply the soundness gates, and assemble the result.

"""
    FunnelData

The (possibly partial) certified funnel of a chain, forward-ordered.

- `ellipsoids` — funnel ellipsoids. On SUCCESS they cover states `1..K+1`. On a
  BACKWARD failure at `failed_k` they cover states `failed_k+1 .. K+1` — i.e.
  `ellipsoids[1]` is the funnel at state `failed_k + 1` (`failed_k` doubles as
  the index offset; `bidirectional_certify!` and `prefix_replan_certify!` rely
  on this). On a FORWARD failure at `failed_k` they cover states `1..failed_k`;
- `kappas` — the matching controllers (`|ellipsoids| = |kappas| + 1` whenever
  at least one step certified);
- `reason` — why the chain could not even start (no terminal/entry ellipsoid,
  endpoint gate failure), or `nothing`.
"""
struct FunnelData{TE, TK}
    ellipsoids::Vector{TE}
    kappas::Vector{TK}
    reason::Union{Nothing, String}
end

FunnelData(ellipsoids, kappas) = FunnelData(ellipsoids, kappas, nothing)

_empty_funnel(reason) = FunnelData(
    LazySets.Ellipsoid{Float64, Vector{Float64}, Matrix{Float64}}[],
    MS.AffineMap[],
    reason,
)

"""
    CertificationResult

Outcome of a certification chain.

- `success` — every step certified *and* every enabled gate passed (including
  terminal containment when checked);
- `failed_k` — first failing step index, or `nothing`. Backward chains fail AT
  step `failed_k` (states `failed_k+1..K+1` stay certified — see
  [`FunnelData`](@ref)); forward chains fail at step `failed_k` with states
  `1..failed_k` certified; `failed_k == K + 1` is a complete chain whose
  terminal gate failed;
- `steps` — forward-ordered [`StepRecord`](@ref)s;
- `controller` — the per-step affine controllers `κ_1..κ_K` on success, else `nothing`;
- `terminal_contained` — whether the terminal ellipsoid lies inside the problem's
  target set (`nothing` when the check is disabled);
- `initial_coverage` — max of `(v − c₁)ᵀP₁(v − c₁)` over the initial set's vertices
  (≤ 1 means the entry funnel covers the initial set; `nothing` if unavailable);
- `state_domain_checked` — whether the reach-avoid gate could run on this domain type;
- `lmi_data` — the [`FunnelData`](@ref).
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

# The ONE assembler both directions use: takes forward-ordered steps and
# ellipsoids (backward chains reverse before calling).
function _assemble_result(
    success,
    failed_k,
    t0,
    steps_forward,
    ellipsoids_forward,
    terminal_contained,
    initial_cov,
    domain_checked;
    reason = nothing,
)
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
        FunnelData(ellipsoids_forward, kappas, reason),
    )
end

function _failed_before_start(failed_k, t0, terminal_contained, reason)
    return CertificationResult(
        false,
        failed_k,
        Float64(time() - t0),
        StepRecord[],
        nothing,
        terminal_contained,
        nothing,
        false,
        _empty_funnel(reason),
    )
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

function run_chain!(ctx::ChainContext)
    t0 = time()
    opts = ctx.options

    E_next, terminal_reason = _terminal_ellipsoid(ctx)
    if E_next !== nothing && terminal_reason === nothing
        # The terminal enters no StepRecord — gate it here or never.
        terminal_reason = endpoint_gate(
            E_next,
            ctx.problem.system.X,
            opts.r_min,
            opts.check_state_domain,
            "terminal",
        )
        terminal_reason === nothing || (E_next = nothing)
    end
    E_next === nothing && return _failed_before_start(ctx.K, t0, false, terminal_reason)

    terminal_contained =
        opts.check_terminal ? terminal_containment(E_next, ctx.problem.target_set) : nothing

    domain_checked =
        opts.check_state_domain && _state_domain_supported(ctx.problem.system.X)

    steps = StepRecord[]
    ellipsoids = [E_next]

    for k in ctx.K:-1:1
        rec = backward_step!(ctx, k, E_next)
        rec = apply_gates(rec, ctx)

        # Two-step rescue: when the one-step transition dies, retry as a single
        # two-step transition into E_{k+2} through the already-synthesized
        # κ_{k+1} — no intermediate containment (see two_step.jl). The funnel
        # invariant then holds across the pair as a composition: from E_k the
        # state is guaranteed in E_{k+2} after two steps (E_{k+1} membership is
        # NOT claimed), and the step-indexed controller replay is unaffected.
        if rec.status != :ok && opts.two_step_rescue
            rec2 = _two_step_rescue!(ctx, k, steps, ellipsoids)
            if rec2 !== nothing
                rec2 = apply_gates(rec2, ctx)
                rec2.status == :ok && (rec = rec2)
            end
        end

        push!(steps, rec)

        if rec.status != :ok
            return _assemble_result(
                false,
                k,
                t0,
                reverse(steps),
                reverse(ellipsoids),
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
        reverse(steps),
        reverse(ellipsoids),
        terminal_contained,
        initial_cov,
        domain_checked,
    )
end
