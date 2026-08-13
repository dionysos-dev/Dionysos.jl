# The shared kernel skeleton. A solver status is not a certificate:
# `ALMOST_OPTIMAL` / `NEARLY_FEASIBLE_POINT` solutions may violate the PSD
# blocks, and then the S-procedure guarantee simply does not hold. Every kernel
# therefore re-validates its constraints numerically at the returned solution.
#
# The skeleton makes that validation structurally drift-proof: `_pose_psd!`
# takes a block *builder* and its argument tuple, poses the JuMP constraint,
# and records a closure that rebuilds the SAME builder on the value-mapped SAME
# arguments — posing and validating cannot disagree. Scalar side-conditions and
# multiplier signs are recorded the same way (`_check!`, `_require_nonneg!`).
#
# `_VALIDATION_TOL` absorbs eigensolver roundoff only — the posed constraints
# demand ⪰ ε·I with ε = `_PSD_STRICTNESS`, so genuinely feasible solutions
# clear 0 with margin.

const _VALIDATION_TOL = 1e-9

# Strictness margin of every posed PSD block. It must DOMINATE the SDP
# solver's feasibility tolerance (Clarabel default 1e-8): the solver returns
# solutions violating posed constraints by up to that tolerance, so an active
# block lands at ε ± tol — with ε ≫ tol it stays PSD and validation only has
# to absorb eigensolver roundoff. At ε ≤ tol the margin is platform roundoff
# luck (measured: same problem validated on Windows, rejected on Linux).
const _PSD_STRICTNESS = 1e-7

# Solution values of numbers, JuMP objects, and (possibly mixed) arrays thereof.
_solved(x::Number) = x
_solved(x::AbstractArray{<:Number}) = x
_solved(x::AbstractArray) = _solved.(x)
_solved(x) = value(x)

# One transition-synthesis SDP under construction: the JuMP model, the
# strictness margin `ε` of its PSD constraints, and the recorded rebuild
# closures / scalar checks that re-validate the certificate at the solution.
# (Comment, not docstring: underscore-internals are excluded from the docs'
# @autodocs blocks, and a docstring outside the manual fails `checkdocs`.)
struct _KernelModel
    model::Model
    ε::Float64
    rebuilds::Vector{Function}
    checks::Vector{Function}
end

_KernelModel(sdp_solver, ε::Float64) =
    _KernelModel(Model(sdp_solver), ε, Function[], Function[])

# Pose `builder(args...) ⪰ ε·I` and record the numeric rebuild of the same call.
function _pose_psd!(km::_KernelModel, builder, args...; ε = km.ε)
    blk = builder(args...)
    @constraint(km.model, blk >= _eye(size(blk, 1)) * ε, PSDCone())
    push!(km.rebuilds, () -> builder(map(_solved, args)...))
    return nothing
end

# Record a scalar validation check (evaluated at the solution).
function _check!(chk::Function, km::_KernelModel)
    push!(km.checks, chk)
    return nothing
end

# S-procedure multipliers must come out nonnegative for the certificate to hold.
function _require_nonneg!(km::_KernelModel, vars...)
    for v in vars
        _check!(() -> all(_solved(v) .>= -_VALIDATION_TOL), km)
    end
    return nothing
end

function _solver_accepted(model::Model)
    term = termination_status(model)
    pstat = primal_status(model)
    return term in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL) &&
           pstat in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
end

# The full a-posteriori validation: every recorded scalar check, then every
# recorded block rebuilt at the solution and required PSD.
function _validated(km::_KernelModel, label::String)
    for chk in km.checks
        if !chk()
            @debug "$label: solver-accepted solution failed a scalar check"
            return false
        end
    end
    for rebuild in km.rebuilds
        M = rebuild()
        margin = LA.eigmin(LA.Symmetric(Matrix{Float64}(M)))
        if margin < -_VALIDATION_TOL
            @debug "$label: solver-accepted solution failed block validation" margin
            return false
        end
    end
    return true
end

# ------------------------------------------------------------
# Remainder models. The linearization error over the box
# [−Lip[1:nx], Lip[1:nx]]·scale enters the reach blocks either as its exact 2ⁿ
# corners (`:vertices` — tightest, exponential block count), as one scalar norm
# ball (`:ball` — constant block count, every axis pays the full radius), or as
# the box's John ellipsoid √n·diag(Lip)·B (`:john_ball` — per-axis radii at
# :ball's block count; measured decisive on the double pendulum, whose position
# rows are ~10× below its velocity rows).
#
# SQUARED-RADIUS CONVENTION: the box scale the kernels feed in is δx + δu where
# δx / δu bound the SQUARED state/input deviations (see `_source_radius_block`
# and `_input_proximity_block`) — the remainder `Lip·(δx + δu)` is a Hessian
# bound, quadratic in the deviation radii.
# ------------------------------------------------------------

struct _RemainderSpec
    ball_like::Bool
    vertices::Union{Nothing, Vector{Vector{Float64}}}
    ρnorm::Float64
    Rsq::Matrix{Float64}
end

# Vertices of the Lipschitz box.
function _lipschitz_vertices(Lip, nx)
    r = collect(Lip[1:nx])
    X = LazySets.Hyperrectangle(zeros(eltype(r), nx), r)
    return collect.(LazySets.vertices_list(X))
end

function _remainder_spec(remainder_model::Symbol, Lip, nx)
    remainder_model in (:vertices, :ball, :john_ball) || error(
        "remainder_model must be :vertices, :ball, or :john_ball, " *
        "got $remainder_model.",
    )
    if remainder_model === :vertices
        return _RemainderSpec(false, _lipschitz_vertices(Lip, nx), NaN, _eye(nx))
    end
    Lipx = collect(Float64, Lip[1:nx])
    ρnorm = remainder_model === :john_ball ? sqrt(nx) : LA.norm(Lipx)
    Rsq = remainder_model === :john_ball ? LA.diagm(Lipx .^ 2) : _eye(nx)
    return _RemainderSpec(true, nothing, ρnorm, Rsq)
end

# Pose the reach blocks for every (remainder direction, noise vertex) pair.
# `base_aux` is the affine part without remainder and noise; `scale` the (fixed
# or variable) squared-deviation scale multiplying the Lipschitz box.
function _pose_reach!(
    km::_KernelModel,
    spec::_RemainderSpec,
    beta,
    muball,
    shape,
    At,
    base_aux,
    νs,
    Q2,
    scale,
    nx,
)
    if spec.ball_like
        for (j, ν) in enumerate(νs)
            aux0 = @expression(km.model, base_aux + hcat(Vector(ν)))
            ρ = @expression(km.model, spec.ρnorm * scale)
            _pose_psd!(
                km,
                _reach_ball_block,
                beta[1, j],
                shape,
                At,
                aux0,
                Q2,
                ρ,
                muball[j],
                nx,
                spec.Rsq,
            )
        end
    else
        for (i, μ) in enumerate(spec.vertices), (j, ν) in enumerate(νs)
            aux = @expression(km.model, base_aux + hcat(μ) * scale + hcat(Vector(ν)))
            _pose_psd!(km, _reach_block, beta[i, j], shape, At, aux, Q2, nx)
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Common constraint groups
# ------------------------------------------------------------

function _pose_inputs!(km::_KernelModel, tau, U, shape, Kmat, ell, nx)
    for i in eachindex(U)
        _pose_psd!(km, _input_block, tau[i], shape, U[i] * Kmat, U[i] * ell, nx)
    end
    return nothing
end

function _pose_cost!(km::_KernelModel, γ, J, shape, G, d, Λ, nx)
    _pose_psd!(km, _cost_block, γ, J, shape, G, d, Λ, nx)
    return nothing
end

function _pose_proximity!(km::_KernelModel, ϕ, shape, F, ell, u_ref, δu, nx, nu)
    _pose_psd!(km, _input_proximity_block, ϕ, shape, F, ell - hcat(u_ref), δu, nx, nu)
    return nothing
end

# Squared-deviation caps: `maxδx`/`maxδu` are RADII at the public boundary, the
# variables bound radii² — hence the squares.
function _pose_deviation_caps!(km::_KernelModel, δx, maxδx, δu, maxδu)
    @constraint(km.model, δx <= maxδx^2)
    @constraint(km.model, δu <= maxδu^2)
    _check!(() -> _solved(δx) <= maxδx^2 + _VALIDATION_TOL, km)
    _check!(() -> _solved(δu) <= maxδu^2 + _VALIDATION_TOL, km)
    return nothing
end

# Per-axis cap on the synthesized source: Q₁ = L·Lᵀ, so Q₁[i,i] = ‖L[i,:]‖² and
# the slab containment |xᵢ − c₁ᵢ| ≤ dᵢ over the ellipsoid is the SOC row
# ‖L[i,:]‖ ≤ dᵢ.
function _pose_source_cap!(km::_KernelModel, source_cap, L, nx)
    source_cap === nothing && return nothing
    @assert length(source_cap) == nx "source_cap must have one entry per state."
    for i in 1:nx
        @constraint(km.model, vcat(source_cap[i], L[i, :]) in SecondOrderCone())
    end
    _check!(km) do
        return all(
            sum(abs2, _solved(L)[i, :]) <= source_cap[i]^2 + _VALIDATION_TOL for i in 1:nx
        )
    end
    return nothing
end

# Size term of the source-synthesis objective (`min λ·J − (1−λ)·size`).
# `:maximin` maximizes the smallest semi-axis (eig(L) are the semi-axes since
# Q = L·Lᵀ): a pancake ellipsoid — the observed collapse mode of long chains —
# is worthless under it, and it needs no exotic cone. `:logdet` is true volume;
# `:trace` its stable proxy.
function _pose_source_size_objective!(km::_KernelModel, objective, λ, L, J, nx)
    if objective === :maximin && λ < 1.0
        r = @variable(km.model, lower_bound = 0.0)
        # The floor is ACTIVE at the optimum (r is maximized onto it), so it
        # must carry the full strictness margin like every validated block —
        # posed at ε = 0 its rebuilt eigmin is pure solver noise and validation
        # becomes a coin flip. The reported floor only understates by ε.
        _pose_psd!(km, _size_floor_block, L, r, nx)
        @objective(km.model, Min, λ * J - (1.0 - λ) * r)
    elseif objective === :logdet && λ < 1.0
        t = @variable(km.model)
        L_tri = [L[i, j] for j in 1:nx for i in j:nx]
        @constraint(km.model, vcat(t, 1.0, L_tri) in MOI.LogDetConeTriangle(nx))
        @objective(km.model, Min, λ * J - (1.0 - λ) * t)
    elseif objective in (:trace, :logdet, :maximin)
        @objective(km.model, Min, λ * J - (1.0 - λ) * sum(L[i, i] for i in 1:nx))
    else
        error("objective must be :maximin, :logdet, or :trace, got $objective.")
    end
    return nothing
end

# Backward-style epilogue: L·Lᵀ is directly the shape matrix Q₁ of the
# synthesized source, and κ = [F·L⁻¹ ℓ].
function _source_solution(L, F, ell, label::String)
    Lval = _solved(L)
    kappa = [_solved(F) / Lval _solved(ell)]
    if !all(isfinite, kappa)
        @debug "$label: near-singular source shape, K = F·L⁻¹ blew up"
        return nothing, nothing
    end
    return Lval * transpose(Lval), kappa
end
