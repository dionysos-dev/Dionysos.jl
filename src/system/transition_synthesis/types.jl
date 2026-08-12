# Result type, the affine-system union, and the set → LMI-data normalization
# every kernel shares.

# n×n dense identity, used across the LMI builders both as identity blocks inside
# matrix literals and as the `ε·_eye` PSD regularization term.
_eye(n) = LA.diagm(ones(n))

"""
    AffineSys

Union of the affine discrete-time system types the transition-synthesis kernels
accept: `(Noisy)ConstrainedAffineControlDiscreteSystem` and
`ConstrainedAffineControlMap`.
"""
AffineSys = Union{
    HybridSystems.NoisyConstrainedAffineControlDiscreteSystem,
    HybridSystems.ConstrainedAffineControlDiscreteSystem,
    HybridSystems.HybridSystems.ConstrainedAffineControlMap,
}

"""
    TransitionResult

Outcome of one transition-synthesis SDP (`solve_transition`,
`solve_transition_backward`, `solve_transition_backward_2step`,
`solve_transition_forward`).

- `feasible::Bool` — whether a certified controller was found;
- `controller` — the affine controller as a `MathematicalSystems.AffineMap`
  `x ↦ Kx + b` (with `b = ℓ − K·c₁`), or `nothing` if infeasible;
- `cost` — upper bound on the worst-case transition cost, or `nothing`;
- `source` — the synthesized source ellipsoid (backward modes only), or `nothing`;
- `target` — the synthesized target ellipsoid (forward mode only), or `nothing`.
"""
struct TransitionResult
    feasible::Bool
    controller::Union{Nothing, MS.AffineMap}
    cost::Union{Nothing, Float64}
    source::Union{Nothing, LazySets.Ellipsoid}
    target::Union{Nothing, LazySets.Ellipsoid}
end

TransitionResult(feasible, controller, cost, source) =
    TransitionResult(feasible, controller, cost, source, nothing)

_infeasible_transition() = TransitionResult(false, nothing, nothing, nothing, nothing)

"""
    as_controller(result::TransitionResult) -> Union{Nothing, AffineController}

The synthesized affine controller as a simulatable [`AffineController`](@ref),
or `nothing` if the transition was infeasible.
"""
as_controller(result::TransitionResult) =
    result.controller === nothing ? nothing : AffineController(result.controller)

# ------------------------------------------------------------
# Argument normalization: sets → LMI data
# ------------------------------------------------------------

"""
    format_input_set(U) -> Vector{<:AbstractMatrix}

Convert the input set `U` (`LazySets.AbstractHyperrectangle`,
`LazySets.Ellipsoid`, or a `LazySets.IntersectionArray` of those) into the
list of matrices `Uᵢ` encoding the input constraints `|Uᵢ·u| ≤ 1` used by the
transition-synthesis LMIs.
"""
function format_input_set(rec::LazySets.AbstractHyperrectangle)
    n = LazySets.dim(rec)
    Uaux = LA.diagm(1:n)
    U = [(Uaux .== i) ./ LazySets.high(rec, i) for i in 1:n]
    return U
end

function format_input_set(elli::LazySets.Ellipsoid)
    return [sqrt(UT.get_quadratic_form(elli))]
end

function format_input_set(iset::LazySets.IntersectionArray)
    return reduce(vcat, [format_input_set(set) for set in iset.array])
end

"""
    format_noise_set(rec::LazySets.AbstractHyperrectangle) -> Matrix

Vertices of the noise polytope `rec` as an `n × 2ⁿ` matrix (one vertex per
column), the format consumed by the transition-synthesis LMIs.
"""
function format_noise_set(rec::LazySets.AbstractHyperrectangle)
    # `vertices_list` repeats vertices on degenerate (zero-radius) axes; each
    # duplicate would add an identical PSD block to every transition SDP.
    return reduce(hcat, unique(LazySets.vertices_list(rec)))
end

_input_matrices(U::AbstractVector) = U
_input_matrices(U) = format_input_set(U)

_noise_vertex_matrix(W::AbstractMatrix) = W
_noise_vertex_matrix(W::LazySets.AbstractHyperrectangle) = format_noise_set(W)

# Noise vertices in state space: systems carrying a noise matrix D map the
# noise-space vertices through it; systems without D take them as given.
_state_noise(::Any, W) = _noise_vertex_matrix(W)
_state_noise(affsys::HybridSystems.NoisyConstrainedAffineControlDiscreteSystem, W) =
    affsys.D * _noise_vertex_matrix(W)

# Every reach block quantifies over the noise vertices — with none, the SDP
# would have NO reachability constraint left and would "certify" any transition.
function _check_noise_vertices(Wcols)
    size(Wcols, 2) >= 1 || error(
        "W must provide at least one disturbance vertex (use a single zero " *
        "vertex for deterministic dynamics) — with none, the reachability " *
        "constraints would vanish and the SDP would certify nothing.",
    )
    return Wcols
end

# The LMIs bound ‖Λ·[x; u; 1]‖² = [x;u;1]'·Λ'Λ·[x;u;1]: they consume a *factor* Λ
# of the PSD cost matrix S = Λ'Λ. Cost conversion is centralized in
# `UT._cost_matrix` (the PSD matrix); the factor is taken here — a symmetric
# square root, so for idempotent S (identity, projections) the factor is S itself
# and the historical golden values are preserved.
function _cost_factor(cost)
    S = LA.Symmetric(Matrix{Float64}(UT._cost_matrix(cost)))
    E = LA.eigen(S)
    return LA.Diagonal(sqrt.(clamp.(E.values, 0.0, Inf))) * transpose(E.vectors)
end

# kappa = [K ℓ] → (K, ℓ)
function _split_controller(kappa)
    nu = size(kappa, 1)
    return kappa[:, 1:(end - 1)], SVector{nu}(kappa[:, end])
end
