# Ellipsoid-to-ellipsoid transition synthesis for affine systems: given
# x⁺ = Ax + Bu + c (+ Dw), decide whether an affine controller
# u(x) = K(x − c₁) + ℓ drives a source ellipsoid into a target ellipsoid under
# input constraints and polytopic noise, while minimizing an upper bound J on
# the worst-case transition cost ‖Λ·[x; u; 1]‖². Each entry point assembles one
# SDP from shared S-procedure blocks; see Corollary 1 of
# https://arxiv.org/pdf/2204.00315.pdf.
#
# Entry points: `solve_transition` (both shapes fixed),
# `solve_transition_backward` (source shape is a decision variable),
# `stabilizing_feedback` (Lyapunov pair for an affine system).

import LazySets
using JuMP

# n×n dense identity, used across the LMI builders both as identity blocks inside
# matrix literals and as the `ε·eye` PSD regularization term.
eye(n) = LA.diagm(ones(n))

AffineSys = Union{
    HybridSystems.NoisyConstrainedAffineControlDiscreteSystem,
    HybridSystems.ConstrainedAffineControlDiscreteSystem,
    HybridSystems.HybridSystems.ConstrainedAffineControlMap,
}

"""
    TransitionResult

Outcome of one transition-synthesis SDP (`solve_transition` /
`solve_transition_backward`).

- `feasible::Bool` — whether a certified controller was found;
- `controller` — the affine controller as a `MathematicalSystems.AffineMap`
  `x ↦ Kx + b` (with `b = ℓ − K·c₁`), or `nothing` if infeasible;
- `cost` — upper bound on the worst-case transition cost, or `nothing`;
- `source` — the synthesized source ellipsoid (backward mode only), or `nothing`;
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
    result = Any[]
    for set in iset.array
        Base.append!(result, format_input_set(set))
    end
    return result
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

# ------------------------------------------------------------
# Shared S-procedure PSD blocks. `shape` is the source-set shape matrix
# (P₁ when the source is fixed, I when its square root is the variable).
#
# Each block has a *builder* returning the matrix — used both to pose the JuMP
# constraint and, with the returned numbers plugged back in, to validate the
# certificate independently of the solver (`_validate_blocks`). Keeping one
# builder per block is what guarantees the validator checks exactly what the
# SDP promised.
# ------------------------------------------------------------

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

_reach_constraint!(model, β, shape, At, aux, P2inv, n, ε) = @constraint(
    model,
    _reach_block(β, shape, At, aux, P2inv, n) >= eye(2 * n + 1) * ε,
    PSDCone()
)

# Input feasibility: |Uᵢ·u(x)| ≤ 1 on the source set, for one constraint
# matrix (UK = Uᵢ·K, Uℓ = Uᵢ·ℓ).
function _input_block(τ, shape, UK, Uℓ, n)
    n_ui = size(UK, 1)
    z = zeros(n, 1)
    return [
        τ*shape z transpose(UK)
        transpose(z) 1-τ transpose(Uℓ)
        UK Uℓ eye(n_ui)
    ]
end

_input_constraint!(model, τ, shape, UK, Uℓ, n, ε) = @constraint(
    model,
    _input_block(τ, shape, UK, Uℓ, n) >= eye(n + size(UK, 1) + 1) * ε,
    PSDCone()
)

# Cost bound: ‖Λ·[x; u(x); 1]‖² ≤ J on the source set, with
# [x; u(x); 1]' = ξ'·G + d for ξ in the source set.
function _cost_block(γ, J, shape, G, d, Λ, n)
    nS = size(Λ, 1)
    z = zeros(n, 1)
    return [
        γ*shape z G*transpose(Λ)
        transpose(z) J-γ d*transpose(Λ)
        Λ*transpose(G) Λ*transpose(d) eye(nS)
    ]
end

_cost_constraint!(model, γ, J, shape, G, d, Λ, n, ε) = @constraint(
    model,
    _cost_block(γ, J, shape, G, d, Λ, n) >= eye(n + size(Λ, 1) + 1) * ε,
    PSDCone()
)

# ------------------------------------------------------------
# A-posteriori certificate validation. A solver status is not a certificate:
# `ALMOST_OPTIMAL` / `NEARLY_FEASIBLE_POINT` solutions may violate the PSD
# blocks, and then the S-procedure guarantee simply does not hold. Every kernel
# therefore rebuilds its blocks numerically at the returned solution and
# requires them PSD. `_VALIDATION_TOL` absorbs eigensolver roundoff only — the
# posed constraints demanded ⪰ ε·I with ε ≥ 1e-10, so genuinely feasible
# solutions clear 0 with margin.
# ------------------------------------------------------------

const _VALIDATION_TOL = 1e-9

function _validate_blocks(blocks)
    for M in blocks
        LA.eigmin(LA.Symmetric(Matrix{Float64}(M))) >= -_VALIDATION_TOL || return false
    end
    return true
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
function _reach_ball_block(β, shape, At, aux0, P2inv, ρ, μ, n)
    z = zeros(n, 1)
    return [
        β*shape z transpose(At) z
        transpose(z) 1-β transpose(aux0) ρ
        At aux0 P2inv-μ*eye(n) z
        transpose(z) ρ transpose(z) μ
    ]
end

# Input proximity: ‖u(x) − u_ref‖² ≤ δu on the source set {ξ : ξ'·shape·ξ ≤ 1}
# (`ellmu` = ℓ − u_ref; `shape` = I when the source square root is the variable,
# P₁ when the source is fixed).
function _input_proximity_block(ϕ, shape, F, ellmu, δu, nx, nu)
    z = zeros(nx, 1)
    return [
        ϕ*shape z transpose(F)
        transpose(z) δu-ϕ transpose(ellmu)
        F ellmu eye(nu)
    ]
end

# Source-radius bound: L·L' ⪯ δx·I.
function _source_radius_block(L, δx, nx)
    return [
        eye(nx) L
        transpose(L) δx*eye(nx)
    ]
end

# ------------------------------------------------------------
# Fixed-shape kernel: both ellipsoids given, decision variables (K, ℓ).
# ------------------------------------------------------------

# The target enters the reachability blocks as its inverse quadratic form,
# which is exactly the LazySets shape matrix Q2 — no inversion needed.
function _fixed_shape_kernel(A, B, c, Wcols, U, Λ, c1, P1, c2, Q2, sdp_solver)
    n = length(c)
    Nw = size(Wcols, 2)
    Nu = length(U)
    m = size(U[1], 2)
    ε = 1e-10

    model = Model(sdp_solver)
    @variable(model, K[i = 1:m, j = 1:n])
    @variable(model, ell[i = 1:m, j = 1:1])
    @variable(model, beta[i = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)

    @expressions(model, begin
        At, A + B * K
        ct, c + B * ell
    end)

    z = zeros(n, 1)

    for i in 1:Nw
        w = Wcols[:, i]
        aux = A * hcat(c1) + hcat(ct) - hcat(c2) + hcat(w)
        _reach_constraint!(model, beta[i], P1, At, aux, Q2, n, ε)
    end

    for i in 1:Nu
        _input_constraint!(model, tau[i], P1, U[i] * K, U[i] * ell, n, ε)
    end

    G = [LA.I transpose(K) z]
    d = [transpose(c1) transpose(ell) 1]
    _cost_constraint!(model, γ, J, P1, G, d, Λ, n, ε)

    @objective(model, Min, J)

    optimize!(model)

    term = termination_status(model)
    pstat = primal_status(model)
    if !(
        term in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL) &&
        pstat in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    )
        return false, nothing, nothing
    end

    Kv = value.(K)
    ellv = value.(ell)
    βv = value.(beta)
    τv = value.(tau)

    Atv = A + B * Kv
    ctv = c + B * ellv
    blocks = Any[]
    for i in 1:Nw
        aux = A * hcat(c1) + hcat(ctv) - hcat(c2) + hcat(Wcols[:, i])
        push!(blocks, _reach_block(βv[i], P1, Atv, aux, Q2, n))
    end
    for i in 1:Nu
        push!(blocks, _input_block(τv[i], P1, U[i] * Kv, U[i] * ellv, n))
    end
    Gv = [LA.I transpose(Kv) z]
    dv = [transpose(c1) transpose(ellv) 1]
    push!(blocks, _cost_block(value(γ), value(J), P1, Gv, dv, Λ, n))

    validated =
        all(βv .>= -_VALIDATION_TOL) &&
        all(τv .>= -_VALIDATION_TOL) &&
        _validate_blocks(blocks)
    if !validated
        @debug "solve_transition: solver-accepted solution failed validation" term pstat
        return false, nothing, nothing
    end

    return true, [Kv ellv], value(J)
end

"""
    solve_transition(affsys, source, target, U, W, cost, sdp_solver) -> TransitionResult

Synthesize an affine controller `u(x) = K(x − c₁) + ℓ` certifying that every
state of the `source` ellipsoid reaches the `target` ellipsoid in one step of
the affine system `affsys` (`x⁺ = Ax + Bu + c (+ Dw)`), for every disturbance
vertex, while minimizing an upper bound on the worst-case transition cost
`‖cost · [x; u; 1]‖²`.

# Arguments
- `affsys`: an affine system (`(Noisy)ConstrainedAffineControlDiscreteSystem`
  or `ConstrainedAffineControlMap`);
- `source`, `target`: `LazySets.Ellipsoid`s;
- `U`: input constraints — a set (`LazySets.AbstractHyperrectangle`, `LazySets.Ellipsoid`,
  `LazySets.IntersectionArray`) or a preformatted list (`format_input_set`);
- `W`: disturbance vertices — one per column, or a box (`LazySets.AbstractHyperrectangle`); mapped
  through the system's noise matrix `D` when the system has one;
- `cost`: the PSD cost matrix `S` of size `(nx+nu+1) × (nx+nu+1)` bounding
  `[x; u; 1]ᵀ·S·[x; u; 1]`, or a `UT.QuadraticStateControlFunction`; the factor
  the LMIs consume is taken internally (`_cost_factor`);
- `sdp_solver`: a JuMP-compatible SDP optimizer (Clarabel, Mosek, …).

Every returned certificate is re-validated numerically at the solver's solution
(`_validate_blocks`) — a transition reported feasible is PSD-certified
independently of the solver's termination status.
"""
function solve_transition(
    affsys::AffineSys,
    source::LazySets.Ellipsoid,
    target::LazySets.Ellipsoid,
    U,
    W,
    cost,
    sdp_solver,
)
    c1 = LazySets.center(source)
    feasible, kappa, J = _fixed_shape_kernel(
        affsys.A,
        affsys.B,
        affsys.c,
        _state_noise(affsys, W),
        _input_matrices(U),
        _cost_factor(cost),
        c1,
        UT.get_quadratic_form(source),
        LazySets.center(target),
        LazySets.shape_matrix(target),
        sdp_solver,
    )
    feasible || return _infeasible_transition()
    K, ℓ = _split_controller(kappa)
    return TransitionResult(true, MS.AffineMap(K, ℓ - K * c1), J, nothing)
end

# ------------------------------------------------------------
# Free-source kernel: target given, the source shape (via its square root L)
# is a decision variable; the source center c₁ and a reference input are data.
# ------------------------------------------------------------

# Vertices of the Lipschitz box: the linearization-error term of the
# reachability blocks scales the box [−Lip[1:nx], Lip[1:nx]] by (δx + δu).
function _lipschitz_vertices(Lip, nx)
    r = collect(Lip[1:nx])
    X = LazySets.Hyperrectangle(zeros(eltype(r), nx), r)
    return LazySets.vertices_list(X)
end

function _free_source_kernel(
    A,
    B,
    c,
    Wcols,
    U,
    Λ,
    c1,
    u_ref,
    c2,
    Q2,
    Lip,
    sdp_solver;
    maxδx,
    maxδu,
    λ,
    objective::Symbol,
    remainder_model::Symbol,
    source_cap = nothing,
)
    nx = length(c)
    nu = size(U[1], 2)
    μ = remainder_model === :ball ? nothing : _lipschitz_vertices(Lip, nx)
    ρnorm = LA.norm(collect(Float64, Lip[1:nx]))
    ν = [Wcols[:, i] for i in 1:size(Wcols, 2)]
    Nx = remainder_model === :ball ? 1 : length(μ)
    Nw = length(ν)
    Nu = length(U)
    ε = 1e-8

    model = Model(sdp_solver)
    @variable(model, L[i = 1:nx, j = 1:nx], PSD)
    @variable(model, F[i = 1:nu, j = 1:nx])

    # Per-axis cap on the source: Q₁ = L·Lᵀ, so Q₁[i,i] = ‖L[i,:]‖² and the slab
    # containment |xᵢ − c₁ᵢ| ≤ dᵢ over the ellipsoid is the SOC row ‖L[i,:]‖ ≤ dᵢ.
    if source_cap !== nothing
        @assert length(source_cap) == nx "source_cap must have one entry per state."
        for i in 1:nx
            @constraint(model, vcat(source_cap[i], L[i, :]) in SecondOrderCone())
        end
    end
    @variable(model, ell[i = 1:nu, j = 1:1])
    @variable(model, beta[i = 1:Nx, j = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, δx >= 0)
    @variable(model, δu >= 0)
    @variable(model, ϕ >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)

    @expressions(model, begin
        At, A * L + B * F
        ct, c + B * ell
    end)

    z = zeros(nx, 1)

    local muball
    if remainder_model === :ball
        @variable(model, muball[j = 1:Nw] >= 0)
        for j in 1:Nw
            aux0 =
                @expression(model, A * hcat(c1) + hcat(ct) - hcat(c2) + hcat(Vector(ν[j])))
            ρ = @expression(model, ρnorm * (δx + δu))
            @constraint(
                model,
                _reach_ball_block(beta[1, j], eye(nx), At, aux0, Q2, ρ, muball[j], nx) >=
                eye(2 * nx + 2) * ε,
                PSDCone()
            )
        end
    else
        for i in 1:Nx, j in 1:Nw
            aux = @expression(
                model,
                A * hcat(c1) + hcat(ct) - hcat(c2) +
                hcat(Vector(μ[i])) * (δx + δu) +
                hcat(Vector(ν[j]))
            )
            _reach_constraint!(model, beta[i, j], eye(nx), At, aux, Q2, nx, ε)
        end
    end

    for i in 1:Nu
        _input_constraint!(model, tau[i], eye(nx), U[i] * F, U[i] * ell, nx, ε)
    end

    G = [transpose(L) transpose(F) z]
    d = [transpose(c1) transpose(ell) 1]
    _cost_constraint!(model, γ, J, eye(nx), G, d, Λ, nx, ε)

    # Input proximity: ‖u(x) − u_ref‖² ≤ δu on the source set.
    u = hcat(u_ref)
    @constraint(
        model,
        _input_proximity_block(ϕ, eye(nx), F, ell - u, δu, nx, nu) >= eye(nx + nu + 1) * ε,
        PSDCone()
    )

    # Source-radius bound: L·L' ⪯ δx·I.
    @constraint(model, _source_radius_block(L, δx, nx) >= eye(nx * 2) * ε, PSDCone())

    @constraint(model, δx <= maxδx^2)
    @constraint(model, δu <= maxδu^2)

    # Volume term of the objective. `:maximin` maximizes the smallest semi-axis
    # (eig(L) are the semi-axes since Q = L·Lᵀ): a pancake ellipsoid — the observed
    # collapse mode of long chains — is worthless under it, and it needs no exotic
    # cone. `:logdet` is true volume; `:trace` its stable proxy.
    r = nothing
    if objective === :maximin && λ < 1.0
        r = @variable(model, lower_bound = 0.0)
        @constraint(model, L >= r * eye(nx), PSDCone())
        @objective(model, Min, λ * J - (1.0 - λ) * r)
    elseif objective === :logdet && λ < 1.0
        @variable(model, t)

        # Lower-triangular entries of symmetric PSD matrix L
        L_tri = [L[i, j] for j in 1:nx for i in j:nx]

        @constraint(model, vcat(t, 1.0, L_tri) in MOI.LogDetConeTriangle(nx),)

        @objective(model, Min, λ * J - (1.0 - λ) * t)
    elseif objective in (:trace, :logdet, :maximin)
        @objective(model, Min, λ * J - (1.0 - λ) * sum(L[i, i] for i in 1:nx),)
    else
        error("objective must be :maximin, :logdet, or :trace, got $objective.")
    end

    optimize!(model)

    term = termination_status(model)
    pstat = primal_status(model)

    if !(
        term in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL) &&
        pstat in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    )
        @debug "solve_transition_backward: infeasible SDP" term pstat raw_status(model)
        return false, nothing, nothing, nothing
    end

    Lval = value.(L)
    Fval = value.(F)
    ellv = value.(ell)
    βv = value.(beta)
    τv = value.(tau)
    δxv = value(δx)
    δuv = value(δu)
    ϕv = value(ϕ)

    Atv = A * Lval + B * Fval
    ctv = c + B * ellv
    blocks = Any[]
    if remainder_model === :ball
        mbv = value.(muball)
        all(mbv .>= -_VALIDATION_TOL) || begin
            @debug "solve_transition_backward: negative ball multiplier"
            return false, nothing, nothing, nothing
        end
        for j in 1:Nw
            aux0 = A * hcat(c1) + hcat(ctv) - hcat(c2) + hcat(Vector(ν[j]))
            ρv = ρnorm * (δxv + δuv)
            push!(
                blocks,
                _reach_ball_block(βv[1, j], eye(nx), Atv, aux0, Q2, ρv, mbv[j], nx),
            )
        end
    else
        for i in 1:Nx, j in 1:Nw
            aux =
                A * hcat(c1) + hcat(ctv) - hcat(c2) +
                hcat(Vector(μ[i])) * (δxv + δuv) +
                hcat(Vector(ν[j]))
            push!(blocks, _reach_block(βv[i, j], eye(nx), Atv, aux, Q2, nx))
        end
    end
    for i in 1:Nu
        push!(blocks, _input_block(τv[i], eye(nx), U[i] * Fval, U[i] * ellv, nx))
    end
    Gv = [transpose(Lval) transpose(Fval) z]
    dv = [transpose(c1) transpose(ellv) 1]
    push!(blocks, _cost_block(value(γ), value(J), eye(nx), Gv, dv, Λ, nx))
    push!(
        blocks,
        _input_proximity_block(ϕv, eye(nx), Fval, ellv - hcat(u_ref), δuv, nx, nu),
    )
    push!(blocks, _source_radius_block(Lval, δxv, nx))
    r === nothing || push!(blocks, Lval - value(r) * eye(nx))

    validated =
        all(βv .>= -_VALIDATION_TOL) &&
        all(τv .>= -_VALIDATION_TOL) &&
        δxv <= maxδx^2 + _VALIDATION_TOL &&
        δuv <= maxδu^2 + _VALIDATION_TOL &&
        (
            source_cap === nothing ||
            all(sum(abs2, Lval[i, :]) <= source_cap[i]^2 + _VALIDATION_TOL for i in 1:nx)
        ) &&
        _validate_blocks(blocks)
    if !validated
        @debug "solve_transition_backward: solver-accepted solution failed validation" term pstat
        return false, nothing, nothing, nothing
    end

    # L·Lᵀ is directly the shape matrix Q of the synthesized source.
    Q1 = Lval * transpose(Lval)
    kappa = [Fval / Lval ellv]
    if !all(isfinite, kappa)
        @debug "solve_transition_backward: near-singular source shape, K = F·L⁻¹ blew up"
        return false, nothing, nothing, nothing
    end
    return true, Q1, kappa, value(J)
end

"""
    solve_transition_backward(affsys, target, source_center, u_ref, U, W, cost,
                              lipschitz, sdp_solver;
                              maxδx = 100.0, maxδu = 20.0, λ = 0.01,
                              use_log_det = true) -> TransitionResult

Synthesize an affine controller *and the largest source ellipsoid* centered at
`source_center` that the controller certifiably drives into the `target`
ellipsoid in one step of `affsys`, robustly to the disturbance vertices `W` and
to the linearization error bounded by `lipschitz` (the Lipschitz radii of
`[x; u; 1]`, scaled by the synthesized state/input deviations `δx + δu`).

The `objective` trades the transition-cost bound against the source size:
`:maximin` maximizes the smallest semi-axis (`min λ·J − (1−λ)·r` with `L ⪰ r·I` —
collapse-proof and cone-free, the certifier's default), `:logdet` the true volume,
`:trace` its stable proxy. (`use_log_det = true/false` is the deprecated spelling of
`:logdet`/`:trace` and wins when passed.)
`u_ref` is the reference input the controller stays `δu`-close to; `maxδx` /
`maxδu` cap the synthesized deviations. `source_cap` (a per-state vector `d`)
additionally confines the source to the axis-aligned slab `|xᵢ − c₁ᵢ| ≤ dᵢ` —
one SOC row per state on the shape factor — so a chain can keep its funnels
inside a state domain by construction. Other arguments as in
[`solve_transition`](@ref). Returns a [`TransitionResult`](@ref) whose `source`
is the synthesized ellipsoid.
"""
function solve_transition_backward(
    affsys::AffineSys,
    target::LazySets.Ellipsoid,
    source_center,
    u_ref,
    U,
    W,
    cost,
    lipschitz,
    sdp_solver;
    maxδx = 100.0,
    maxδu = 20.0,
    λ = 0.01,
    objective::Symbol = :logdet,
    remainder_model::Symbol = :vertices,
    source_cap = nothing,
    use_log_det = nothing,
)
    # `use_log_det` is the pre-maximin spelling; it wins when passed explicitly.
    use_log_det !== nothing && (objective = use_log_det ? :logdet : :trace)
    feasible, Q1, kappa, J = _free_source_kernel(
        affsys.A,
        affsys.B,
        affsys.c,
        _state_noise(affsys, W),
        _input_matrices(U),
        _cost_factor(cost),
        source_center,
        u_ref,
        LazySets.center(target),
        LazySets.shape_matrix(target),
        lipschitz,
        sdp_solver;
        maxδx = maxδx,
        maxδu = maxδu,
        λ = λ,
        objective = objective,
        remainder_model = remainder_model,
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

# ------------------------------------------------------------
# Free-target kernel (forward direction): the source ellipsoid is data — so the
# linearization deviation δx is a *constant* and the remainder is known before
# solving — and the target enters the reach blocks only as its shape matrix Q₂,
# where it is linear. A free target is the naturally convex direction of ellipsoid
# calculus: no L-congruence, no volume cone, no adaptive boxes (plan.md §4.4).
# ------------------------------------------------------------

function _free_target_kernel(
    A,
    B,
    c,
    Wcols,
    U,
    Λ,
    c1,
    P1,
    δx_const,
    u_ref,
    c2,
    Qhat,          # target shape matrix for the α-mode; `nothing` = free-shape mode
    Lip,
    sdp_solver;
    maxδu,
    λ,
    q_min,
    q_max,
    remainder_model::Symbol,
)
    nx = length(c)
    nu = size(U[1], 2)
    μ = remainder_model === :ball ? nothing : _lipschitz_vertices(Lip, nx)
    ρnorm = LA.norm(collect(Float64, Lip[1:nx]))
    ν = [Wcols[:, i] for i in 1:size(Wcols, 2)]
    Nx = remainder_model === :ball ? 1 : length(μ)
    Nw = length(ν)
    Nu = length(U)
    ε = 1e-8

    model = Model(sdp_solver)
    @variable(model, K[i = 1:nu, j = 1:nx])
    @variable(model, ell[i = 1:nu, j = 1:1])
    @variable(model, beta[i = 1:Nx, j = 1:Nw] >= 0)
    @variable(model, tau[i = 1:Nu] >= 0)
    @variable(model, δu >= 0)
    @variable(model, ϕ >= 0)
    @variable(model, γ >= 0)
    @variable(model, J >= 0)

    local α, Q2
    if Qhat !== nothing
        @variable(model, α >= 0)
        Q2 = α .* Qhat
    else
        @variable(model, Q2[i = 1:nx, j = 1:nx], PSD)
        # Conditioning sandwich: thin targets are fine for the certificate but the
        # next step inverts Q₂ into its source form.
        @constraint(model, Q2 >= q_min * eye(nx), PSDCone())
        @constraint(model, q_max * eye(nx) >= Q2, PSDCone())
    end

    @expressions(model, begin
        At, A + B * K
        ct, c + B * ell
    end)

    z = zeros(nx, 1)

    local muball
    if remainder_model === :ball
        @variable(model, muball[j = 1:Nw] >= 0)
        for j in 1:Nw
            aux0 =
                @expression(model, A * hcat(c1) + hcat(ct) - hcat(c2) + hcat(Vector(ν[j])))
            ρ = @expression(model, ρnorm * (δx_const + δu))
            @constraint(
                model,
                _reach_ball_block(beta[1, j], P1, At, aux0, Q2, ρ, muball[j], nx) >=
                eye(2 * nx + 2) * ε,
                PSDCone()
            )
        end
    else
        for i in 1:Nx, j in 1:Nw
            aux = @expression(
                model,
                A * hcat(c1) + hcat(ct) - hcat(c2) +
                hcat(Vector(μ[i])) * (δx_const + δu) +
                hcat(Vector(ν[j]))
            )
            _reach_constraint!(model, beta[i, j], P1, At, aux, Q2, nx, ε)
        end
    end

    for i in 1:Nu
        _input_constraint!(model, tau[i], P1, U[i] * K, U[i] * ell, nx, ε)
    end

    G = [LA.I transpose(K) z]
    d = [transpose(c1) transpose(ell) 1]
    _cost_constraint!(model, γ, J, P1, G, d, Λ, nx, ε)

    # Input proximity feeds the u-part of the remainder: ‖u(x) − u_ref‖² ≤ δu.
    u = hcat(u_ref)
    @constraint(
        model,
        _input_proximity_block(ϕ, P1, K, ell - u, δu, nx, nu) >= eye(nx + nu + 1) * ε,
        PSDCone()
    )

    @constraint(model, δu <= maxδu^2)

    size_term = Qhat !== nothing ? α : sum(Q2[i, i] for i in 1:nx)
    @objective(model, Min, λ * J + (1.0 - λ) * size_term)

    optimize!(model)

    term = termination_status(model)
    pstat = primal_status(model)
    if !(
        term in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL) &&
        pstat in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    )
        @debug "solve_transition_forward: infeasible SDP" term pstat raw_status(model)
        return false, nothing, nothing, nothing
    end

    Kv = value.(K)
    ellv = value.(ell)
    βv = value.(beta)
    τv = value.(tau)
    δuv = value(δu)
    ϕv = value(ϕ)
    Q2v = Qhat !== nothing ? value(α) .* Qhat : value.(Q2)

    Atv = A + B * Kv
    ctv = c + B * ellv
    blocks = Any[]
    if remainder_model === :ball
        mbv = value.(muball)
        all(mbv .>= -_VALIDATION_TOL) || begin
            @debug "solve_transition_forward: negative ball multiplier"
            return false, nothing, nothing, nothing
        end
        for j in 1:Nw
            aux0 = A * hcat(c1) + hcat(ctv) - hcat(c2) + hcat(Vector(ν[j]))
            ρv = ρnorm * (δx_const + δuv)
            push!(blocks, _reach_ball_block(βv[1, j], P1, Atv, aux0, Q2v, ρv, mbv[j], nx))
        end
    else
        for i in 1:Nx, j in 1:Nw
            aux =
                A * hcat(c1) + hcat(ctv) - hcat(c2) +
                hcat(Vector(μ[i])) * (δx_const + δuv) +
                hcat(Vector(ν[j]))
            push!(blocks, _reach_block(βv[i, j], P1, Atv, aux, Q2v, nx))
        end
    end
    for i in 1:Nu
        push!(blocks, _input_block(τv[i], P1, U[i] * Kv, U[i] * ellv, nx))
    end
    Gv = [LA.I transpose(Kv) z]
    dv = [transpose(c1) transpose(ellv) 1]
    push!(blocks, _cost_block(value(γ), value(J), P1, Gv, dv, Λ, nx))
    push!(blocks, _input_proximity_block(ϕv, P1, Kv, ellv - hcat(u_ref), δuv, nx, nu))

    validated =
        all(βv .>= -_VALIDATION_TOL) &&
        all(τv .>= -_VALIDATION_TOL) &&
        δuv <= maxδu^2 + _VALIDATION_TOL &&
        _validate_blocks(blocks)
    if !validated
        @debug "solve_transition_forward: solver-accepted solution failed validation" term pstat
        return false, nothing, nothing, nothing
    end

    return true, Q2v, [Kv ellv], value(J)
end

"""
    solve_transition_forward(affsys, source, target_center, u_ref, U, W, cost,
                             lipschitz, sdp_solver;
                             target_shape = nothing, maxδu = 20.0, λ = 0.01,
                             q_min = 1e-9, q_max = 1e9) -> TransitionResult

Synthesize an affine controller *and the smallest target ellipsoid* centered at
`target_center` that the controller certifiably reaches from the **given** `source`
ellipsoid in one step of `affsys`, robustly to the disturbance vertices and to the
linearization error bounded by `lipschitz` — whose state deviation is the *known*
`λ_max(Q₁)` of the source, so the remainder is fixed data (no adaptive boxes).

Target modes: with `target_shape` (a LazySets shape matrix `Q̂`) only the scale `α`
is free (`min λ·J + (1−λ)·α` — `α` is the per-step contraction number); without it
the shape `Q₂` itself is a decision variable (`min λ·J + (1−λ)·tr(Q₂)`) inside the
conditioning sandwich `q_min·I ⪯ Q₂ ⪯ q_max·I`. Both are single convex SDPs — the
target enters the reach blocks linearly. Certificates are validated at the returned
solution as in [`solve_transition`](@ref). Returns a [`TransitionResult`](@ref)
whose `target` is the synthesized ellipsoid.
"""
function solve_transition_forward(
    affsys::AffineSys,
    source::LazySets.Ellipsoid,
    target_center,
    u_ref,
    U,
    W,
    cost,
    lipschitz,
    sdp_solver;
    target_shape = nothing,
    maxδu = 20.0,
    λ = 0.01,
    q_min = 1e-9,
    q_max = 1e9,
    remainder_model::Symbol = :vertices,
)
    c1 = LazySets.center(source)
    P1 = Matrix{Float64}(UT.get_quadratic_form(source))
    δx_const = LA.eigmax(LA.Symmetric(Matrix{Float64}(LazySets.shape_matrix(source))))

    feasible, Q2, kappa, J = _free_target_kernel(
        affsys.A,
        affsys.B,
        affsys.c,
        _state_noise(affsys, W),
        _input_matrices(U),
        _cost_factor(cost),
        c1,
        P1,
        δx_const,
        u_ref,
        target_center,
        target_shape,
        lipschitz,
        sdp_solver;
        maxδu = maxδu,
        λ = λ,
        q_min = q_min,
        q_max = q_max,
        remainder_model = remainder_model,
    )
    feasible || return _infeasible_transition()
    K, ℓ = _split_controller(kappa)
    controller = MS.AffineMap(K, ℓ - K * c1)
    target = LazySets.Ellipsoid(
        collect(float.(target_center)),
        UT._symmetrize(Matrix{Float64}(Q2));
        check_posdef = false,
    )
    return TransitionResult(true, controller, J, nothing, target)
end

# ------------------------------------------------------------
# Stabilizing feedback
# ------------------------------------------------------------

"""
    stabilizing_feedback(subsys, sdp_solver) -> (feasible, K, P, γ)

For a stabilizable affine system, find the state-feedback gain `K` and the
matrix `P` satisfying the discrete-time Lyapunov inequality
`(A+BK)'P(A+BK) − P ≺ 0`, minimizing the condition number of `P` (`γ` is the
attained smallest eigenvalue of `P⁻¹`).
"""
function stabilizing_feedback(
    subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem,
    sdp_solver,
)
    A = subsys.A
    B = subsys.B
    n = size(A, 1)
    m = size(B, 2)

    model = Model(sdp_solver)
    @variable(model, L[i = 1:m, j = 1:n])
    @variable(model, S[i = 1:n, j = 1:n], PSD)
    @variable(model, gamma)

    t(x) = transpose(x)

    @constraint(model, [
        S t(A * S + B * L)
        A * S+B * L S
    ] >= 1e-10 * eye(2n), PSDCone())
    @constraint(model, eye(n) >= S, PSDCone())
    @constraint(model, S >= gamma * eye(n), PSDCone())

    @objective(model, Max, gamma)

    optimize!(model)

    P = inv(value.(S))
    K = value.(L) * P
    gamma = value(gamma)
    feasible = solution_summary(model).termination_status == MOI.OPTIMAL
    return feasible, K, P, gamma
end
