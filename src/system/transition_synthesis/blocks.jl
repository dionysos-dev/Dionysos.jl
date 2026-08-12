# Shared S-procedure PSD blocks. `shape` is the source-set shape matrix
# (P₁ when the source is fixed, I when its square root is the variable).
#
# Each block is a *pure builder* returning the matrix. The kernel skeleton
# (`_pose_psd!`) uses one and the same argument tuple to pose the JuMP
# constraint AND to rebuild the numeric block for the a-posteriori validation —
# so the validator provably checks exactly what the SDP promised, arguments
# included.

# Reachability: the closed loop maps {ξ: ξ'·shape·ξ ≤ 1} into the target
# ellipsoid, for one disturbance vertex (`aux` collects the affine part).
function _reach_block(β, shape, At, aux, P2inv, n)
    z = zeros(n, 1)
    return [
        β*shape z transpose(At)
        transpose(z) 1-β transpose(aux)
        At aux P2inv
    ]
end

# Ball-remainder reach block (plan.md §4.4-★2): the linearization error enters as
# one norm-bounded uncertainty ‖e‖₂ ≤ ρ (the circumscribing ball of the Lipschitz
# box — mildly conservative) instead of 2ⁿ box vertices. By Petersen's lemma the
# robust condition M₀ + q·eᵀ·Pᵀ + P·e·qᵀ ⪰ 0 ∀‖e‖ ≤ ρ (with P selecting the
# target block-row and q the middle row) is, Schur-completed to stay linear in ρ
# and the multiplier μ:
#
#     [ M₀ − μ·P·Pᵀ   ρ·q ]
#     [ ρ·qᵀ           μ  ]  ⪰ 0
#
# — one extra multiplier and one border row per noise vertex, constant in n.
# A shaped uncertainty e ∈ D·{‖ẽ‖ ≤ ρ} (the `:john_ball` model) substitutes
# e = D·ẽ, which turns P into P·D and hence μ·P·Pᵀ into μ·P·D²·Pᵀ — the
# target-block deduction becomes μ·Rsq with Rsq = D², still linear in μ.
function _reach_ball_block(β, shape, At, aux0, P2inv, ρ, μ, n, Rsq = _eye(n))
    z = zeros(n, 1)
    return [
        β*shape z transpose(At) z
        transpose(z) 1-β transpose(aux0) ρ
        At aux0 P2inv-μ*Rsq z
        transpose(z) ρ transpose(z) μ
    ]
end

# Input feasibility: |Uᵢ·u(x)| ≤ 1 on the source set, for one constraint
# matrix (UK = Uᵢ·K, Uℓ = Uᵢ·ℓ).
function _input_block(τ, shape, UK, Uℓ, n)
    n_ui = size(UK, 1)
    z = zeros(n, 1)
    return [
        τ*shape z transpose(UK)
        transpose(z) 1-τ transpose(Uℓ)
        UK Uℓ _eye(n_ui)
    ]
end

# Cost bound: ‖Λ·[x; u(x); 1]‖² ≤ J on the source set, with
# [x; u(x); 1]' = ξ'·G + d for ξ in the source set.
function _cost_block(γ, J, shape, G, d, Λ, n)
    nS = size(Λ, 1)
    z = zeros(n, 1)
    return [
        γ*shape z G*transpose(Λ)
        transpose(z) J-γ d*transpose(Λ)
        Λ*transpose(G) Λ*transpose(d) _eye(nS)
    ]
end

# Input proximity: ‖u(x) − u_ref‖² ≤ δu on the source set {ξ : ξ'·shape·ξ ≤ 1}
# (`ellmu` = ℓ − u_ref; `shape` = I when the source square root is the variable,
# P₁ when the source is fixed). NOTE the squared-radius convention: δu bounds the
# SQUARED input deviation, hence the `maxδu^2` caps at the kernel level.
function _input_proximity_block(ϕ, shape, F, ellmu, δu, nx, nu)
    z = zeros(nx, 1)
    return [
        ϕ*shape z transpose(F)
        transpose(z) δu-ϕ transpose(ellmu)
        F ellmu _eye(nu)
    ]
end

# Source-radius bound: L·L' ⪯ δx·I — δx bounds the SQUARED source semi-axes
# (same squared-radius convention as δu above).
function _source_radius_block(L, δx, nx)
    return [
        _eye(nx) L
        transpose(L) δx*_eye(nx)
    ]
end

# Maximin floor: L ⪰ r·I, so maximizing r maximizes the smallest semi-axis.
_size_floor_block(L, r, nx) = L - r * _eye(nx)
