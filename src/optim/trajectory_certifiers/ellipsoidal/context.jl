# The chain context: everything `run_chain!` needs, built once per certification —
# problem, trajectory, affine-approximation provider, transition cost, SDP backend,
# options — plus the optional per-step state scaling and the transition-solve wrapper.

struct ChainContext{P, TR, AP, TX, TU, TS, TB, OPTS}
    problem::P
    traj::TR
    affine_provider::AP
    xs::TX
    us::TU
    K::Int
    S::TS
    backend::TB
    options::OPTS
end

function _identity_transition_cost(nx::Int, nu::Int)
    return [
        Matrix{Float64}(LA.I, nx, nx) zeros(nx, nu) zeros(nx, 1)
        zeros(nu, nx) Matrix{Float64}(LA.I, nu, nu) zeros(nu, 1)
        zeros(1, nx + nu) 1.0
    ]
end

function _transition_cost_matrix(cost, nx, nu)
    S =
        cost === nothing ? _identity_transition_cost(nx, nu) :
        Matrix{Float64}(UT._cost_matrix(cost))
    @assert size(S) == (nx + nu + 1, nx + nu + 1)
    return S
end

function build_context(
    problem::PR.ProblemType,
    traj::ST.Trajectory,
    affine_provider,
    backend,
    options::ChainOptions,
)
    xs = collect(ST.states(traj))
    us = collect(ST.inputs(traj))

    K = length(us)
    @assert length(xs) == K + 1 "Expected length(xs) == length(us) + 1."

    nx = length(xs[1])
    nu = length(us[1])
    S = _transition_cost_matrix(options.transition_cost, nx, nu)

    return ChainContext(problem, traj, affine_provider, xs, us, K, S, backend, options)
end

# ------------------------------------------------------------
# Scaling helpers (per-step change of coordinates x = xk + T·z)
# ------------------------------------------------------------

# Coordinate change z = T⁻¹(x − xnext). In quadratic form Pz = TᵀPT, which in
# shape-matrix form is Qz = T⁻¹·Q·T⁻ᵀ — no inversion of Q needed.
function _scale_target_ellipsoid(E, xnext, scaling)
    scaling === nothing && return E

    Tinv = inv(Matrix{Float64}(scaling))
    Qz = UT._symmetrize(Tinv * Matrix(LazySets.shape_matrix(E)) * Tinv')
    cz = collect(Tinv * (collect(LazySets.center(E)) - collect(xnext)))

    return LazySets.Ellipsoid(cz, Qz; check_posdef = false)
end

# Inverse change x = T·z + center: Qx = T·Qz·Tᵀ.
function _unscale_source_ellipsoid(Ez, center, scaling)
    scaling === nothing && return Ez

    T = Matrix{Float64}(scaling)
    Qx = UT._symmetrize(T * Matrix(LazySets.shape_matrix(Ez)) * T')

    return LazySets.Ellipsoid(collect(float.(center)), Qx; check_posdef = false)
end

function _scaled_affine_system(affsys, xk, xnext, scaling)
    scaling === nothing && return affsys

    Tstate = Matrix{Float64}(scaling)
    Tinv = inv(Tstate)

    A = Matrix(Tinv * affsys.A * Tstate)
    B = Matrix(Tinv * affsys.B)
    c = vec(Tinv * (affsys.A * xk + affsys.c - xnext))
    Dnoise = Matrix(Tinv * affsys.D)

    return MS.NoisyConstrainedAffineControlDiscreteSystem(
        A,
        B,
        c,
        Dnoise,
        affsys.X,
        affsys.U,
        affsys.W,
    )
end

# CAVEAT (plan.md §4.2-B): this maps the error *vector* to z-coordinates but leaves
# the (δx + δu) multiplier in z-units, while the Hessian bound was taken in
# x-coordinates — sound only when σ_max(scaling) ≤ 1 (contractive scaling). A
# magnifying scaling under-estimates the remainder by up to σ_max². Until the
# provider computes bounds on the scaled dynamics (plan P3), refuse magnifying
# scalings outright.
function _scaled_lipschitz(L, nx::Int, scaling)
    scaling === nothing && return L

    T = Matrix{Float64}(scaling)
    σmax = LA.opnorm(T)
    σmax <= 1.0 + 1e-12 || error(
        "state_scaling has σ_max = $σmax > 1: the scaled Lipschitz remainder would " *
        "be under-estimated by up to σ_max² (unsound). Normalize the dynamics " *
        "globally instead, or use a contractive scaling.",
    )

    Ls = collect(Float64, L)
    Tinv = inv(T)
    Ls[1:nx] .= abs.(Tinv) * Ls[1:nx]

    return Ls
end

# Cost matrices are PSD ([x; u; 1]ᵀ·S·[x; u; 1]); under x = xk + T·z the transform is
# the congruence Tᵀ·S·T (the kernels take the factor themselves).
function _scaled_transition_cost(S, xk, nu::Int, scaling)
    scaling === nothing && return S

    nx = length(xk)
    T = zeros(nx + nu + 1, nx + nu + 1)

    T[1:nx, 1:nx] .= Matrix{Float64}(scaling)
    T[1:nx, end] .= xk
    T[(nx + 1):(nx + nu), (nx + 1):(nx + nu)] .= Matrix{Float64}(LA.I, nu, nu)
    T[end, end] = 1.0

    return transpose(T) * S * T
end

# ------------------------------------------------------------
# Transition solve
# ------------------------------------------------------------

function _transition_backward(
    affsys,
    E_next,
    xk,
    xnext,
    uk,
    Uformat,
    Wformat,
    S,
    L,
    backend;
    maxδx,
    maxδu,
    λ,
    state_scaling,
    objective = :maximin,
    remainder_model = :vertices,
    source_cap = nothing,
)
    if state_scaling === nothing
        result = ST.solve_transition_backward(
            affsys,
            E_next,
            xk,
            uk,
            Uformat,
            Wformat,
            S,
            L,
            backend;
            maxδx = maxδx,
            maxδu = maxδu,
            λ = λ,
            objective = objective,
            remainder_model = remainder_model,
            source_cap = source_cap,
        )
        return result.source, result.controller, result.cost
    end

    nx = length(xk)

    affsys_z = _scaled_affine_system(affsys, xk, xnext, state_scaling)
    E_next_z = _scale_target_ellipsoid(E_next, xnext, state_scaling)
    Lz = _scaled_lipschitz(L, nx, state_scaling)
    Sz = _scaled_transition_cost(S, xk, length(uk), state_scaling)

    # The kernel's slab cap is axis-aligned; only a diagonal per-step scaling maps
    # an axis-aligned slab in x to one in z.
    source_cap_z = nothing
    if source_cap !== nothing
        T = Matrix{Float64}(state_scaling)
        LA.isdiag(T) || error(
            "domain_cap with a non-diagonal state_scaling is unsupported: the " *
            "state-domain slab does not stay axis-aligned in scaled coordinates.",
        )
        source_cap_z = collect(Float64, source_cap) ./ LA.diag(T)
    end

    result = ST.solve_transition_backward(
        affsys_z,
        E_next_z,
        zeros(nx),
        uk,
        Uformat,
        Wformat,
        Sz,
        Lz,
        backend;
        maxδx = maxδx,
        maxδu = maxδu,
        λ = λ,
        objective = objective,
        remainder_model = remainder_model,
        source_cap = source_cap_z,
    )

    result.feasible || return nothing, nothing, nothing

    E_prev = _unscale_source_ellipsoid(result.source, xk, state_scaling)

    # source center is the origin in scaled coordinates, so controller.c == ℓ
    Kz = Matrix{Float64}(result.controller.A)
    ell = vec(Float64.(result.controller.c))
    Kx = Kz * inv(Matrix{Float64}(state_scaling))
    cont = MS.AffineMap(Kx, ell - Kx * xk)

    return E_prev, cont, result.cost
end

# Per-step slab cap from the state domain: the distance from the nominal center to
# each domain face, shaved so solver-tolerance solutions still pass the state-domain
# gate. `nothing` when the domain is unsupported or the nominal sits outside it (the
# gate then reports the violation as before).
function _domain_cap(X, xk)
    box = X isa UT.SetMinus ? UT.minus_included(X) : X
    box isa LazySets.AbstractHyperrectangle || return nothing
    cap = min.(LazySets.high(box) .- xk, xk .- LazySets.low(box))
    all(cap .> 0) || return nothing
    return cap .* (1 - 1e-6)
end

function _solve_transition(
    ctx::ChainContext,
    approx,
    E_next,
    xk,
    xnext,
    uk;
    box_cap = nothing,
)
    source_cap =
        ctx.options.domain_cap ? _domain_cap(ctx.problem.system.X, collect(Float64, xk)) :
        nothing
    # Cap the source to the linearization box as well: state-side box consistency
    # holds by construction, so the adaptive search line-searches box scales for
    # the largest certifiable funnel instead of chasing a size-maximizing SDP with
    # ever-bigger boxes (whose Hessian bounds eventually kill the LMI — measured
    # `lmi_infeasible_at_max_box` on the double pendulum's mid-swing).
    if box_cap !== nothing
        shaved = collect(Float64, box_cap) .* (1 - 1e-6)
        source_cap = source_cap === nothing ? shaved : min.(source_cap, shaved)
    end
    return _transition_backward(
        approx.system,
        E_next,
        xk,
        xnext,
        uk,
        approx.Uformat,
        approx.Wformat,
        ctx.S,
        approx.lipschitz,
        ctx.backend;
        maxδx = ctx.options.maxδx,
        maxδu = ctx.options.maxδu,
        λ = ctx.options.λ,
        state_scaling = ctx.options.state_scaling,
        objective = ctx.options.objective,
        remainder_model = ctx.options.remainder_model,
        source_cap = source_cap,
    )
end
