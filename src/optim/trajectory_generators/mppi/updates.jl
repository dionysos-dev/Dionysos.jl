# Update rules over the sampled batch: softmin weighting (:mppi) and elite refit
# (:cem), plus the ESS-adaptive temperature (plan.md §3.1-E1).

# Softmin weights at temperature λ. Returns (weights, ess, entropy).
function _softmin_weights(costs::Vector{Float64}, λ::Float64)
    β = minimum(costs)
    w = exp.(-(costs .- β) ./ λ)
    Z = sum(w)
    if !isfinite(Z) || Z <= 0
        fill!(w, 1 / length(w))
    else
        w ./= Z
    end
    ess = 1 / sum(abs2, w)
    entropy = -sum(x -> x > 0 ? x * log(x) : 0.0, w)
    return w, ess, entropy
end

# Pick λ so the effective sample size hits `ess_target` (in samples): ESS(λ) is
# monotone increasing in λ, so bisect on log λ around the cost spread. This is what
# keeps the weights meaningful regardless of the cost scale — a fixed λ against
# costs of the wrong magnitude degenerates into argmin (weight collapse).
function _ess_lambda(costs::Vector{Float64}, ess_target::Float64)
    scale = max(maximum(costs) - minimum(costs), 1e-12)
    lo, hi = 1e-8 * scale, 1e8 * scale
    for _ in 1:60
        mid = sqrt(lo * hi)
        _, ess, _ = _softmin_weights(costs, mid)
        if ess < ess_target
            lo = mid
        else
            hi = mid
        end
    end
    return sqrt(lo * hi)
end

# :mppi — weighted average of the *effective* noise (what was actually simulated
# after projection; plan.md §3.1-E3), re-projected.
function _mppi_new_controls(u_nom, eff, weights, project_input)
    horizon = length(u_nom)
    return map(1:horizon) do k
        δu = weights[1] * eff[1][k]
        for s in 2:length(weights)
            δu = δu + weights[s] * eff[s][k]
        end
        return project_input(u_nom[k] + δu)
    end
end

# :cem — mean of the elite fraction's applied controls, re-projected.
function _cem_new_controls(u_nom, eff, costs, elite_frac, project_input)
    ns = length(costs)
    nelite = max(2, round(Int, elite_frac * ns))
    elite = partialsortperm(costs, 1:nelite)
    horizon = length(u_nom)
    return map(1:horizon) do k
        acc = u_nom[k] + eff[elite[1]][k]
        for s in elite[2:end]
            acc = acc + u_nom[k] + eff[s][k]
        end
        return project_input(acc / nelite)
    end
end
