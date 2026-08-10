# Noise models. The structured Gaussian type is what unlocks the theory-side
# features (importance-sampling correction, annealing, covariance adaptation); a
# legacy `(rng, u, k) -> ε` closure keeps working as the escape hatch — it simply
# opts out of those features.

"""
    GaussianMPPINoise(σ)

Zero-mean diagonal Gaussian exploration noise with per-dimension standard
deviations `σ` (an `SVector` keeps the sampling loop allocation-free). Enables the
importance-sampling cross term (`α < 1`) and σ-annealing (`anneal`).
"""
struct GaussianMPPINoise{S}
    σ::S
    Σinv_diag::S
end

function GaussianMPPINoise(σ::AbstractVector)
    Σinv = map(s -> s > 0 ? 1 / s^2 : 0.0, σ)
    return GaussianMPPINoise(σ, Σinv)
end

_draw(n::GaussianMPPINoise{<:SVector}, rng, u, k, scale) =
    scale .* n.σ .* randn(rng, typeof(n.σ))
_draw(n::GaussianMPPINoise, rng, u, k, scale) = scale .* n.σ .* randn(rng, length(n.σ))

# Legacy closure: draws are whatever the user's sampler returns; `scale` (annealing)
# does not apply.
_draw(sampler, rng, u, k, scale) = sampler(rng, u, k)

_supports_is_term(::GaussianMPPINoise) = true
_supports_is_term(::Any) = false

# uᵀΣ⁻¹ε for the importance-sampling correction.
_is_dot(n::GaussianMPPINoise, u, ε) = sum(n.Σinv_diag .* u .* ε)
