# Adaptive linearization-box search (backward direction): the SDP synthesizes the
# source ellipsoid, so the box the Hessian bound was taken on is only known to be
# large enough a posteriori. Grow on infeasibility, grow to the required radii on
# inconsistency, then optionally line-search box scales for the largest volume.

_clamp_box_radii(δ, δmin, δmax) = min.(max.(δ, δmin), δmax)
_grow_infeasible_box_radii(δ, δmax, growth) = min.(growth .* δ, δmax)
_grow_to_required_box_radii(required, δmin, δmax, safety) =
    min.(max.(safety .* required, δmin), δmax)

function _evaluate_adaptive_box_candidate(ctx, E_next, k, xk, xnext, uk, δx, δu; atol)
    approx = ST.build_affine_approximation(ctx.affine_provider, xk, uk; δx = δx, δu = δu)

    E_prev, kappa, cost = _solve_transition(ctx, approx, E_next, xk, xnext, uk)

    if E_prev === nothing || kappa === nothing
        return (;
            status = :lmi_infeasible,
            approx,
            E_prev = nothing,
            kappa = nothing,
            cost = nothing,
            required_X_radius = nothing,
            required_U_radius = nothing,
            logvolume = nothing,
        )
    end

    required_X_radius, required_U_radius =
        _required_linearization_box_radii(E_prev, kappa, xk, uk)

    if !_box_contains_required_radii(
        δx,
        δu,
        required_X_radius,
        required_U_radius;
        atol = atol,
    )
        return (;
            status = :inconsistent_box,
            approx,
            E_prev,
            kappa,
            cost,
            required_X_radius,
            required_U_radius,
            logvolume = nothing,
        )
    end

    logvolume = _ellipsoid_logvolume(E_prev)

    if !isfinite(logvolume)
        return (;
            status = :invalid_logvolume,
            approx,
            E_prev,
            kappa,
            cost,
            required_X_radius,
            required_U_radius,
            logvolume = nothing,
        )
    end

    return (;
        status = :ok,
        approx,
        E_prev,
        kappa,
        cost,
        required_X_radius,
        required_U_radius,
        logvolume,
    )
end

function _adaptive_backward_step!(ctx::ChainContext, k::Int, E_next)
    opts = ctx.options.adaptive_boxes
    @assert opts !== nothing "adaptive_boxes cannot be nothing in adaptive mode."

    xk = collect(ctx.xs[k])
    xnext = collect(LazySets.center(E_next))
    uk = collect(ctx.us[k])

    δx = _clamp_box_radii(copy(opts.ΔX_initial), opts.ΔX_min, opts.ΔX_max)
    δu = _clamp_box_radii(copy(opts.ΔU_initial), opts.ΔU_min, opts.ΔU_max)

    last_result = nothing
    last_status = :not_started
    last_iter = 0
    base = nothing

    for iter in 1:opts.max_iters
        result = _evaluate_adaptive_box_candidate(
            ctx,
            E_next,
            k,
            xk,
            xnext,
            uk,
            δx,
            δu;
            atol = opts.atol,
        )

        last_result = result
        last_status = result.status
        last_iter = iter

        if result.status == :lmi_infeasible
            new_δx = _grow_infeasible_box_radii(δx, opts.ΔX_max, opts.growth)
            new_δu = _grow_infeasible_box_radii(δu, opts.ΔU_max, opts.growth)

            if all(new_δx .<= δx .+ opts.atol) && all(new_δu .<= δu .+ opts.atol)
                last_status = :lmi_infeasible_at_max_box
                break
            end

            δx, δu = new_δx, new_δu
            continue
        end

        if result.status == :ok
            base = (; δx = copy(δx), δu = copy(δu), result, iter)

            if opts.keep_first_consistent || opts.objective == :first_consistent
                return StepRecord(
                    k,
                    :ok,
                    Float64(result.cost),
                    result.E_prev,
                    result.kappa,
                    _step_summary(
                        result.approx.lipschitz,
                        δx,
                        δu,
                        result.required_X_radius,
                        result.required_U_radius,
                        iter,
                        :accepted;
                        selected_logvolume = result.logvolume,
                        selected_scale = 1.0,
                        selected_candidate_index = 0,
                        number_of_candidate_boxes = 0,
                        provider_summary = result.approx.summary,
                        _candidate_diagnostics_empty()...,
                    ),
                )
            end

            break
        end

        new_δx = _grow_to_required_box_radii(
            result.required_X_radius,
            opts.ΔX_min,
            opts.ΔX_max,
            opts.safety,
        )
        new_δu = _grow_to_required_box_radii(
            result.required_U_radius,
            opts.ΔU_min,
            opts.ΔU_max,
            opts.safety,
        )

        if all(new_δx .<= δx .+ opts.atol) && all(new_δu .<= δu .+ opts.atol)
            last_status = :inconsistent_at_max_box
            break
        end

        δx, δu = new_δx, new_δu
    end

    if base !== nothing
        diag = (;
            scales = Float64[],
            logvolumes = Union{Nothing, Float64}[],
            statuses = Symbol[],
            Xbar_radii = Vector{Float64}[],
            Ubar_radii = Vector{Float64}[],
        )

        best = base
        best_index = 0
        best_scale = 1.0
        best_logvolume = best.result.logvolume === nothing ? -Inf : best.result.logvolume

        for (idx, scale) in enumerate(opts.search_scales)
            δx_candidate = _clamp_box_radii(scale .* base.δx, opts.ΔX_min, opts.ΔX_max)
            δu_candidate = _clamp_box_radii(scale .* base.δu, opts.ΔU_min, opts.ΔU_max)

            candidate = _evaluate_adaptive_box_candidate(
                ctx,
                E_next,
                k,
                xk,
                xnext,
                uk,
                δx_candidate,
                δu_candidate;
                atol = opts.atol,
            )

            _append_candidate_diagnostic!(
                diag,
                scale,
                δx_candidate,
                δu_candidate,
                candidate,
            )

            if candidate.status == :ok &&
               candidate.logvolume !== nothing &&
               candidate.logvolume > best_logvolume
                best = (;
                    δx = copy(δx_candidate),
                    δu = copy(δu_candidate),
                    result = candidate,
                    iter = base.iter,
                )
                best_index = idx
                best_scale = Float64(scale)
                best_logvolume = candidate.logvolume
            end
        end

        diag_tuple = _candidate_diagnostics_tuple(diag)

        return StepRecord(
            k,
            :ok,
            Float64(best.result.cost),
            best.result.E_prev,
            best.result.kappa,
            _step_summary(
                best.result.approx.lipschitz,
                best.δx,
                best.δu,
                best.result.required_X_radius,
                best.result.required_U_radius,
                best.iter,
                best_index == 0 ? :base_fallback : :max_volume_selected;
                selected_logvolume = best.result.logvolume,
                selected_scale = best_scale,
                selected_candidate_index = best_index,
                number_of_candidate_boxes = length(opts.search_scales),
                provider_summary = best.result.approx.summary,
                diag_tuple...,
            ),
        )
    end

    approx = last_result === nothing ? nothing : last_result.approx

    return StepRecord(
        k,
        :infeasible,
        nothing,
        nothing,
        nothing,
        _step_summary(
            approx === nothing ? nothing : approx.lipschitz,
            δx,
            δu,
            last_result === nothing ? nothing : last_result.required_X_radius,
            last_result === nothing ? nothing : last_result.required_U_radius,
            last_iter,
            last_status;
            provider_summary = approx === nothing ? nothing : approx.summary,
            _candidate_diagnostics_empty()...,
        ),
    )
end
