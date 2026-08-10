# One backward step: certify the transition x̄_k → E_next, producing the step's
# ellipsoid, controller, and diagnostics as a StepRecord. Gates are applied by the
# chain, not here — steps only produce records.

"""
    StepRecord

Outcome of one chain step: index `k`, `status` (`:ok`, `:infeasible`, or a gate
failure), transition-cost bound, certified `ellipsoid`, controller `kappa`, and the
diagnostics `summary` NamedTuple.
"""
struct StepRecord{TE, TK, TS}
    k::Int
    status::Symbol
    cost::Union{Nothing, Float64}
    ellipsoid::TE
    kappa::TK
    summary::TS
end

function StepRecord(
    k::Integer,
    status::Symbol,
    cost,
    ellipsoid,
    kappa,
    summary::TS,
) where {TS}
    return StepRecord(
        Int(k),
        status,
        cost === nothing ? nothing : Float64(cost),
        ellipsoid,
        kappa,
        summary,
    )
end

function _fixed_backward_step!(ctx::ChainContext, k::Int, E_next)
    opts = ctx.options

    @assert !isempty(opts.linearization_δx) "Set options.linearization_δx for fixed mode."
    @assert !isempty(opts.linearization_δu) "Set options.linearization_δu for fixed mode."

    xk = collect(ctx.xs[k])
    xnext = collect(LazySets.center(E_next))
    uk = collect(ctx.us[k])

    δx = collect(Float64, opts.linearization_δx)
    δu = collect(Float64, opts.linearization_δu)

    approx = ST.build_affine_approximation(ctx.affine_provider, xk, uk; δx = δx, δu = δu)

    E_prev, kappa, cost = _solve_transition(ctx, approx, E_next, xk, xnext, uk)

    if E_prev === nothing || kappa === nothing
        return StepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            _step_summary(
                approx.lipschitz,
                δx,
                δu,
                nothing,
                nothing,
                1,
                :fixed_infeasible;
                provider_summary = approx.summary,
                _candidate_diagnostics_empty()...,
            ),
        )
    end

    required_X_radius, required_U_radius =
        _required_linearization_box_radii(E_prev, kappa, xk, uk)

    return StepRecord(
        k,
        :ok,
        Float64(cost),
        E_prev,
        kappa,
        _step_summary(
            approx.lipschitz,
            δx,
            δu,
            required_X_radius,
            required_U_radius,
            1,
            :fixed;
            provider_summary = approx.summary,
            _candidate_diagnostics_empty()...,
        ),
    )
end

function backward_step!(ctx::ChainContext, k::Int, E_next)
    opts = ctx.options.adaptive_boxes
    if opts !== nothing && opts.enabled
        return _adaptive_backward_step!(ctx, k, E_next)
    end

    return _fixed_backward_step!(ctx, k, E_next)
end
