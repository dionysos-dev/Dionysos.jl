export build_symbolic_context, backward_step!, run_backward_chain!

import Dionysos
import IntervalArithmetic as IA
import LinearAlgebra as LA
import MathematicalSystems as MS

const DI = Dionysos
const ST = DI.System
const SY = DI.Symbolic
const UT = DI.Utils

struct AdaptiveLinearizationBoxOptions
    enabled::Bool
    ΔX_initial::Vector{Float64}
    ΔX_min::Vector{Float64}
    ΔX_max::Vector{Float64}
    ΔU_initial::Vector{Float64}
    ΔU_min::Vector{Float64}
    ΔU_max::Vector{Float64}
    growth::Float64
    safety::Float64
    max_iters::Int
    atol::Float64
    verbose::Bool
    search_scales::Vector{Float64}
    objective::Symbol
    keep_first_consistent::Bool
end

struct EllipsoidalBackwardContext{P, C, CFG, SYM, TX, TU, TS, TB}
    problem::P
    candidate::C
    config::CFG
    symbolic::SYM
    xs::TX
    us::TU
    K::Int
    S::TS
    backend::TB
    maxδx::Float64
    maxδu::Float64
    λ::Float64
    terminal_radius::Float64
    state_scaling::Union{Nothing, Vector{Float64}}
    terminal_center::Union{Nothing, Vector{Float64}}
    terminal_shape::Union{Nothing, Matrix{Float64}}
    adaptive_boxes::AdaptiveLinearizationBoxOptions
end

function _identity_transition_cost(nx::Int, nu::Int) # l'utilisateur devrait lui meme fournir la transi cost fonction (pour le moment c'est pratiqque)
    return [
        Matrix{Float64}(LA.I, nx, nx) zeros(nx, nu) zeros(nx, 1)
        zeros(nu, nx) Matrix{Float64}(LA.I, nu, nu) zeros(nu, 1)
        zeros(1, nx + nu) 1.0
    ]
end

function build_symbolic_context(problem, candidate, config, symbolic_builder) # je pourrais cancel ça et tout condenser
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    K = length(us)
    @assert length(xs) == K + 1

    sym = symbolic_builder(problem, candidate, config)

    nx = length(xs[1])
    nu = length(us[1])
    S = _identity_transition_cost(nx, nu)

    opts = config.options
    maxδx = Float64(opts.maxδx)
    maxδu = Float64(opts.maxδu)
    λ = Float64(opts.λ)
    terminal_radius = Float64(opts.rayon_terminal)
    state_scaling = _state_scaling(opts, nx)
    terminal_center, terminal_shape = _terminal_ellipsoid_data(opts, nx)
    adaptive_boxes = _adaptive_linearization_box_options(opts, sym, nx, nu)

    return EllipsoidalBackwardContext(
        problem,
        candidate,
        config,
        sym,
        xs,
        us,
        K,
        S,
        config.backend,
        maxδx,
        maxδu,
        λ,
        terminal_radius,
        state_scaling,
        terminal_center,
        terminal_shape,
        adaptive_boxes,
    )
end

_option(opts, name::Symbol, default) = hasproperty(opts, name) ? getproperty(opts, name) : default

function Base.getproperty(rec::BackwardStepRecord, name::Symbol)
    name === :info && return getfield(rec, :summary)
    return getfield(rec, name)
end

function Base.propertynames(rec::BackwardStepRecord; private::Bool = false)
    names = fieldnames(typeof(rec))
    return private ? (names..., :info) : (names..., :info)
end

function _interval_box_radii(Δ, n::Int, name::Symbol)
    vals = collect(Δ)
    length(vals) == n ||
        throw(ArgumentError("$name must have length $n, got $(length(vals))"))
    radii = [max(abs(Float64(IA.inf(I))), abs(Float64(IA.sup(I)))) for I in vals]
    all(isfinite, radii) ||
        throw(ArgumentError("$name entries must be finite"))
    all(>=(0.0), radii) ||
        throw(ArgumentError("$name entries must be nonnegative"))
    return radii
end

function _vector_option(opts, name::Symbol, default::AbstractVector, n::Int; allow_inf::Bool = false)
    raw = _option(opts, name, default)
    v = Float64.(collect(raw))
    length(v) == n ||
        throw(ArgumentError("$name must have length $n, got $(length(v))"))
    if allow_inf
        all(!isnan, v) ||
            throw(ArgumentError("$name entries must not be NaN"))
    else
        all(isfinite, v) ||
            throw(ArgumentError("$name entries must be finite"))
    end
    all(>=(0.0), v) ||
        throw(ArgumentError("$name entries must be nonnegative"))
    return v
end

function _vector_option(opts, name::Symbol, default::AbstractVector; allow_inf::Bool = false)
    raw = _option(opts, name, default)
    v = Float64.(collect(raw))
    if allow_inf
        all(!isnan, v) ||
            throw(ArgumentError("$name entries must not be NaN"))
    else
        all(isfinite, v) ||
            throw(ArgumentError("$name entries must be finite"))
    end
    all(>=(0.0), v) ||
        throw(ArgumentError("$name entries must be nonnegative"))
    return v
end

function _adaptive_linearization_box_options(opts, sym, nx::Int, nu::Int)
    default_ΔX = _interval_box_radii(sym.ΔX, nx, :ΔX)
    default_ΔU = _interval_box_radii(sym.ΔU, nu, :ΔU)
    enabled = Bool(_option(opts, :adaptive_linearization_boxes, false))

    ΔX_initial = _vector_option(opts, :ΔX_initial, default_ΔX, nx)
    ΔX_min = _vector_option(opts, :ΔX_min, zeros(nx), nx)
    ΔX_max = _vector_option(opts, :ΔX_max, fill(Inf, nx), nx; allow_inf = true)
    ΔU_initial = _vector_option(opts, :ΔU_initial, default_ΔU, nu)
    ΔU_min = _vector_option(opts, :ΔU_min, zeros(nu), nu)
    ΔU_max = _vector_option(opts, :ΔU_max, fill(Inf, nu), nu; allow_inf = true)

    all(ΔX_max .>= ΔX_min) ||
        throw(ArgumentError("ΔX_max must be componentwise >= ΔX_min"))
    all(ΔU_max .>= ΔU_min) ||
        throw(ArgumentError("ΔU_max must be componentwise >= ΔU_min"))

    growth = Float64(_option(opts, :adaptive_box_growth, 2.0))
    safety = Float64(_option(opts, :adaptive_box_safety, 1.05))
    max_iters = Int(_option(opts, :adaptive_box_max_iters, 1))
    atol = Float64(_option(opts, :adaptive_box_atol, 1.0e-8))
    verbose = Bool(_option(opts, :adaptive_box_verbose, false))
    search_scales = _vector_option(
        opts,
        :adaptive_box_search_scales,
        [0.7, 0.85, 1.0, 1.15, 1.35, 1.6, 2.0],
    )
    objective = Symbol(_option(opts, :adaptive_box_objective, :max_volume))
    keep_first_consistent = Bool(_option(opts, :adaptive_box_keep_first_consistent, false))

    isfinite(growth) && growth > 1.0 ||
        throw(ArgumentError("adaptive_box_growth must be finite and > 1"))
    isfinite(safety) && safety >= 1.0 ||
        throw(ArgumentError("adaptive_box_safety must be finite and >= 1"))
    max_iters >= 1 ||
        throw(ArgumentError("adaptive_box_max_iters must be >= 1"))
    isfinite(atol) && atol >= 0.0 ||
        throw(ArgumentError("adaptive_box_atol must be finite and nonnegative"))
    !isempty(search_scales) ||
        throw(ArgumentError("adaptive_box_search_scales must not be empty"))
    all(>(0.0), search_scales) ||
        throw(ArgumentError("adaptive_box_search_scales entries must be strictly positive"))
    objective in (:first_consistent, :max_volume) ||
        throw(ArgumentError("adaptive_box_objective must be :first_consistent or :max_volume"))

    return AdaptiveLinearizationBoxOptions(
        enabled,
        max.(ΔX_initial, ΔX_min),
        ΔX_min,
        ΔX_max,
        max.(ΔU_initial, ΔU_min),
        ΔU_min,
        ΔU_max,
        growth,
        safety,
        max_iters,
        atol,
        verbose,
        search_scales,
        objective,
        keep_first_consistent,
    )
end

function _state_scaling(opts, nx::Int)
    raw = _option(opts, :state_scaling, nothing)
    raw === nothing && return nothing

    s = Float64.(collect(raw))
    length(s) == nx ||
        throw(ArgumentError("state_scaling must have length $nx, got $(length(s))"))
    all(isfinite, s) ||
        throw(ArgumentError("state_scaling entries must be finite"))
    all(>(0.0), s) ||
        throw(ArgumentError("state_scaling entries must be strictly positive"))
    return s
end

function _terminal_ellipsoid_data(opts, nx::Int)
    raw_center = _option(opts, :terminal_center, nothing)
    raw_shape = _option(opts, :terminal_shape, nothing)

    (raw_center === nothing) == (raw_shape === nothing) ||
        throw(ArgumentError("terminal_center and terminal_shape must be provided together"))
    raw_center === nothing && return nothing, nothing

    c = Float64.(collect(raw_center))
    length(c) == nx ||
        throw(ArgumentError("terminal_center must have length $nx, got $(length(c))"))
    all(isfinite, c) ||
        throw(ArgumentError("terminal_center entries must be finite"))

    P = Matrix{Float64}(raw_shape)
    size(P) == (nx, nx) ||
        throw(ArgumentError("terminal_shape must have size ($nx, $nx), got $(size(P))"))
    all(isfinite, P) ||
        throw(ArgumentError("terminal_shape entries must be finite"))
    LA.norm(P - P', Inf) <= 1.0e-9 ||
        throw(ArgumentError("terminal_shape must be symmetric"))
    LA.isposdef(LA.Symmetric(P)) ||
        throw(ArgumentError("terminal_shape must be positive definite"))

    return c, Matrix(LA.Symmetric(P))
end

function _scale_target_ellipsoid(E, xnext, scaling)
    scaling === nothing && return E
    # Convention used here: E(c, P) = {x : (x-c)' P (x-c) <= 1}.
    # With x = xnext + D*z, the normalized ellipsoid center is
    # cz = D^{-1} * (c - xnext), not necessarily zero.
    D = LA.Diagonal(scaling)
    S = LA.Diagonal(1.0 ./ scaling)
    Pz = Matrix(D' * Matrix(E.P) * D)
    cz = collect(S * (collect(E.c) - collect(xnext)))
    return UT.Ellipsoid(Pz, cz)
end

function _unscale_source_ellipsoid(Ez, center, scaling)
    scaling === nothing && return Ez
    # Convention used here: E(c, P) = {x : (x-c)' P (x-c) <= 1}.
    # Since z = D^{-1}(x-c), the physical ellipsoid has Px = D^{-T} * Pz * D^{-1}.
    S = LA.Diagonal(1.0 ./ scaling)
    Px = Matrix(S' * Matrix(Ez.P) * S)
    return UT.Ellipsoid(Px, center)
end

function _scaled_affine_system(affsys, xk, xnext, scaling)
    scaling === nothing && return affsys

    Dstate = LA.Diagonal(scaling)
    Sstate = LA.Diagonal(1.0 ./ scaling)
    A = Matrix(Sstate * affsys.A * Dstate)
    B = Matrix(Sstate * affsys.B)
    c = vec(Sstate * (affsys.A * xk + affsys.c - xnext))
    Dnoise = Matrix(Sstate * affsys.D)
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

function _scaled_lipschitz(L, nx::Int, scaling)
    scaling === nothing && return L
    Ls = collect(Float64, L)
    # The first nx entries are interpreted as output-state remainder radii.
    # They are scaled componentwise because z_next = D^{-1}(x_next - xnext_nominal).
    Ls[1:nx] .= Ls[1:nx] ./ scaling
    return Ls
end

function _scaled_transition_cost(S, xk, nu::Int, scaling)
    scaling === nothing && return S

    # In transition_backward, S is used as a linear factor S * [x; u; 1].
    # With x = xk + D*z, use [x; u; 1] = T * [z; u; 1], hence Sz = S*T.
    nx = length(xk)
    T = zeros(nx + nu + 1, nx + nu + 1)
    T[1:nx, 1:nx] .= LA.Diagonal(scaling)
    T[1:nx, end] .= xk
    T[(nx + 1):(nx + nu), (nx + 1):(nx + nu)] .= Matrix{Float64}(LA.I, nu, nu)
    T[end, end] = 1.0
    return S * T
end

function _transition_backward_with_optional_scaling(
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
)
    state_scaling === nothing && return UT.transition_backward(
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
    )

    nx = length(xk)
    affsys_z = _scaled_affine_system(affsys, xk, xnext, state_scaling)
    E_next_z = _scale_target_ellipsoid(E_next, xnext, state_scaling)
    Lz = _scaled_lipschitz(L, nx, state_scaling)
    Sz = _scaled_transition_cost(S, xk, length(uk), state_scaling)

    Pz, kappa_z, cost = UT.transition_backward(
        affsys_z.A,
        affsys_z.B,
        affsys_z.c,
        affsys_z.D,
        E_next_z.c,
        E_next_z.P,
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
    )

    if Pz === nothing || kappa_z === nothing
        return nothing, nothing, nothing
    end

    E_prev = _unscale_source_ellipsoid(UT.Ellipsoid(Pz, zeros(nx)), xk, state_scaling)
    Kz, ell = UT.get_controller_matrices(kappa_z)
    Kx = Kz * LA.Diagonal(1.0 ./ state_scaling)
    cont = MS.AffineMap(Kx, ell - Kx * xk)
    return E_prev, cont, cost
end

function _centered_interval_box(center::AbstractVector, radius::AbstractVector)
    c = Float64.(collect(center))
    r = Float64.(collect(radius))
    length(c) == length(r) ||
        throw(ArgumentError("center and radius must have the same length"))
    all(isfinite, c) ||
        throw(ArgumentError("box center entries must be finite"))
    all(isfinite, r) ||
        throw(ArgumentError("box radius entries must be finite"))
    all(>=(0.0), r) ||
        throw(ArgumentError("box radius entries must be nonnegative"))
    return IA.IntervalBox(c .- r, c .+ r)
end

function _ellipsoid_axis_radii(E)
    P = Matrix{Float64}(E.P)
    Q = inv(LA.Symmetric(P))
    return sqrt.(max.(0.0, LA.diag(Q)))
end

function _ellipsoid_logvolume(E)
    P = Matrix{Float64}(E.P)
    F = LA.cholesky(LA.Symmetric(P))
    return -sum(log, LA.diag(F.U))
end

function _controller_matrices(kappa, nx::Int)
    if hasproperty(kappa, :A) && hasproperty(kappa, :b)
        return Matrix{Float64}(getproperty(kappa, :A)), vec(Float64.(getproperty(kappa, :b)))
    end
    if hasproperty(kappa, :A) && hasproperty(kappa, :c)
        return Matrix{Float64}(getproperty(kappa, :A)), vec(Float64.(getproperty(kappa, :c)))
    end
    if kappa isa AbstractMatrix
        K = Matrix{Float64}(kappa[:, 1:nx])
        b = vec(Float64.(kappa[:, nx + 1]))
        return K, b
    end
    K, b = UT.get_controller_matrices(kappa)
    return Matrix{Float64}(K), vec(Float64.(b))
end

function _controller_image_axis_radii(kappa, E)
    nx = length(E.c)
    K, b = _controller_matrices(kappa, nx)
    Q = inv(LA.Symmetric(Matrix{Float64}(E.P)))
    Σu = K * Q * K'
    uc = K * collect(Float64, E.c) + b
    ru = sqrt.(max.(0.0, LA.diag(LA.Symmetric(Σu))))
    return uc, ru
end

function _required_linearization_box_radii(E_prev, kappa, xk, uk)
    rx = _ellipsoid_axis_radii(E_prev)
    uc, ru = _controller_image_axis_radii(kappa, E_prev)
    required_δx = abs.(collect(Float64, E_prev.c) .- collect(Float64, xk)) .+ rx
    required_δu = abs.(uc .- collect(Float64, uk)) .+ ru
    return required_δx, required_δu
end

function _box_contains_required_radii(
    δx::AbstractVector,
    δu::AbstractVector,
    required_δx::AbstractVector,
    required_δu::AbstractVector;
    atol::Float64,
)
    return all(required_δx .<= δx .+ atol) && all(required_δu .<= δu .+ atol)
end

_clamp_box_radii(δ::AbstractVector, δmin::AbstractVector, δmax::AbstractVector) =
    min.(max.(δ, δmin), δmax)

_grow_infeasible_box_radii(
    δ::AbstractVector,
    δmax::AbstractVector,
    growth::Float64,
) = min.(growth .* δ, δmax)

_grow_to_required_box_radii(
    required::AbstractVector,
    δmin::AbstractVector,
    δmax::AbstractVector,
    safety::Float64,
) = min.(max.(safety .* required, δmin), δmax)

function _build_affine_approximation(ctx, xk, uk, wk, Xbar, Ubar, Wbar)
    return ST.buildAffineApproximation( #
        ctx.symbolic.fsymbolic,
        ctx.symbolic.x,
        ctx.symbolic.u,
        ctx.symbolic.w,
        xk,
        uk,
        wk,
        Xbar,
        Ubar,
        Wbar,
    )
end

function _solve_backward_transition(ctx, affineSys, E_next, xk, xnext, uk, L)
    return _transition_backward_with_optional_scaling(
        affineSys,
        E_next,
        xk,
        xnext,
        uk,
        ctx.symbolic.Uformat,
        ctx.symbolic.Wformat,
        ctx.S,
        L,
        ctx.backend;
        maxδx = ctx.maxδx,
        maxδu = ctx.maxδu,
        λ = ctx.λ,
        state_scaling = ctx.state_scaling,
    )
end

function _fixed_backward_step!(ctx::EllipsoidalBackwardContext, k::Int, E_next)
    xk = collect(ctx.xs[k])
    xnext = collect(ctx.xs[k + 1])
    uk = collect(ctx.us[k])
    wk = zeros(length(ctx.symbolic.w))

    Xbar = IA.IntervalBox(xk .+ ctx.symbolic.ΔX)
    Ubar = IA.IntervalBox(uk .+ ctx.symbolic.ΔU)
    Wbar = IA.IntervalBox(wk .+ ctx.symbolic.ΔW)

    affineSys, L = _build_affine_approximation(ctx, xk, uk, wk, Xbar, Ubar, Wbar)

    E_prev, kappa, cost =
        _solve_backward_transition(ctx, affineSys, E_next, xk, xnext, uk, L)
    Xbar_radius = _interval_box_radii(ctx.symbolic.ΔX, length(xk), :ΔX)
    Ubar_radius = _interval_box_radii(ctx.symbolic.ΔU, length(uk), :ΔU)

    if E_prev === nothing || kappa === nothing
        println("My transi is impossible")
        return BackwardStepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            (;
                L,
                Xbar_radius,
                Ubar_radius,
                required_X_radius = nothing,
                required_U_radius = nothing,
                adaptive_box_iters = 1,
                adaptive_box_status = :fixed_infeasible,
                selected_logvolume = nothing,
                selected_scale = nothing,
                selected_candidate_index = nothing,
                number_of_candidate_boxes = 0,
                candidate_scales = Float64[],
                candidate_logvolumes = Union{Nothing, Float64}[],
                candidate_statuses = Symbol[],
                candidate_Xbar_radii = Vector{Float64}[],
                candidate_Ubar_radii = Vector{Float64}[],
            ),
        )
    end

    required_X_radius, required_U_radius =
        _required_linearization_box_radii(E_prev, kappa, xk, uk)

    return BackwardStepRecord(
        k,
        :ok,
        Float64(cost),
        E_prev,
        kappa,
        (;
            L,
            Xbar_radius,
            Ubar_radius,
            required_X_radius,
            required_U_radius,
            adaptive_box_iters = 1,
            adaptive_box_status = :fixed,
            selected_logvolume = nothing,
            selected_scale = nothing,
            selected_candidate_index = nothing,
            number_of_candidate_boxes = 0,
            candidate_scales = Float64[],
            candidate_logvolumes = Union{Nothing, Float64}[],
            candidate_statuses = Symbol[],
            candidate_Xbar_radii = Vector{Float64}[],
            candidate_Ubar_radii = Vector{Float64}[],
        ),
    )
end

function _candidate_diagnostics_empty()
    return (;
        candidate_scales = Float64[],
        candidate_logvolumes = Union{Nothing, Float64}[],
        candidate_statuses = Symbol[],
        candidate_Xbar_radii = Vector{Float64}[],
        candidate_Ubar_radii = Vector{Float64}[],
    )
end

function _record_info(
    L,
    δx,
    δu,
    required_X_radius,
    required_U_radius,
    adaptive_box_iters,
    adaptive_box_status;
    selected_logvolume = nothing,
    selected_scale = nothing,
    selected_candidate_index = nothing,
    number_of_candidate_boxes = 0,
    candidate_scales = Float64[],
    candidate_logvolumes = Union{Nothing, Float64}[],
    candidate_statuses = Symbol[],
    candidate_Xbar_radii = Vector{Float64}[],
    candidate_Ubar_radii = Vector{Float64}[],
)
    return (;
        L,
        Xbar_radius = copy(δx),
        Ubar_radius = copy(δu),
        required_X_radius,
        required_U_radius,
        adaptive_box_iters,
        adaptive_box_status,
        selected_logvolume,
        selected_scale,
        selected_candidate_index,
        number_of_candidate_boxes,
        candidate_scales,
        candidate_logvolumes,
        candidate_statuses,
        candidate_Xbar_radii,
        candidate_Ubar_radii,
    )
end

function _evaluate_adaptive_box_candidate(
    ctx,
    E_next,
    xk,
    xnext,
    uk,
    wk,
    Wbar,
    δx,
    δu;
    atol,
)
    Xbar = _centered_interval_box(xk, δx)
    Ubar = _centered_interval_box(uk, δu)
    affineSys, L = _build_affine_approximation(ctx, xk, uk, wk, Xbar, Ubar, Wbar)
    E_prev, kappa, cost =
        _solve_backward_transition(ctx, affineSys, E_next, xk, xnext, uk, L)

    if E_prev === nothing || kappa === nothing
        return (;
            status = :lmi_infeasible,
            L,
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
            L,
            E_prev,
            kappa,
            cost,
            required_X_radius,
            required_U_radius,
            logvolume = nothing,
        )
    end

    logvolume = try
        _ellipsoid_logvolume(E_prev)
    catch
        nothing
    end
    if logvolume === nothing || !isfinite(logvolume)
        return (;
            status = :invalid_logvolume,
            L,
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
        L,
        E_prev,
        kappa,
        cost,
        required_X_radius,
        required_U_radius,
        logvolume,
    )
end

function _append_candidate_diagnostic!(diag, scale, δx, δu, result)
    push!(diag.scales, Float64(scale))
    push!(diag.logvolumes, result.logvolume)
    push!(diag.statuses, result.status)
    push!(diag.Xbar_radii, copy(δx))
    push!(diag.Ubar_radii, copy(δu))
    return diag
end

function _candidate_diagnostics_tuple(diag)
    return (;
        candidate_scales = copy(diag.scales),
        candidate_logvolumes = copy(diag.logvolumes),
        candidate_statuses = copy(diag.statuses),
        candidate_Xbar_radii = copy(diag.Xbar_radii),
        candidate_Ubar_radii = copy(diag.Ubar_radii),
    )
end

function _adaptive_backward_step!(ctx::EllipsoidalBackwardContext, k::Int, E_next)
    opts = ctx.adaptive_boxes
    xk = collect(ctx.xs[k])
    xnext = collect(ctx.xs[k + 1])
    uk = collect(ctx.us[k])
    wk = zeros(length(ctx.symbolic.w))
    Wbar = IA.IntervalBox(wk .+ ctx.symbolic.ΔW)

    δx = _clamp_box_radii(copy(opts.ΔX_initial), opts.ΔX_min, opts.ΔX_max)
    δu = _clamp_box_radii(copy(opts.ΔU_initial), opts.ΔU_min, opts.ΔU_max)

    last_L = nothing
    last_required_X = nothing
    last_required_U = nothing
    last_status = :not_started
    last_iter = 0
    base = nothing

    for iter in 1:opts.max_iters
        result = _evaluate_adaptive_box_candidate(
            ctx,
            E_next,
            xk,
            xnext,
            uk,
            wk,
            Wbar,
            δx,
            δu;
            atol = opts.atol,
        )
        last_iter = iter
        last_L = result.L
        last_status = result.status

        if result.status == :lmi_infeasible
            opts.verbose && println(
                "[adaptive-box] k=",
                k,
                " iter=",
                iter,
                " infeasible δx=",
                δx,
                " δu=",
                δu,
                " L=",
                result.L,
            )
            new_δx = _grow_infeasible_box_radii(δx, opts.ΔX_max, opts.growth)
            new_δu = _grow_infeasible_box_radii(δu, opts.ΔU_max, opts.growth)
            if all(new_δx .<= δx .+ opts.atol) && all(new_δu .<= δu .+ opts.atol)
                last_status = :lmi_infeasible_at_max_box
                break
            end
            δx, δu = new_δx, new_δu
            continue
        end

        last_required_X = result.required_X_radius
        last_required_U = result.required_U_radius

        if result.status == :ok || result.status == :invalid_logvolume
            opts.verbose && println(
                "[adaptive-box] k=",
                k,
                " base accepted iter=",
                iter,
                " δx=",
                δx,
                " δu=",
                δu,
                " required_δx=",
                result.required_X_radius,
                " required_δu=",
                result.required_U_radius,
                " L=",
                result.L,
                " logvolume=",
                result.logvolume,
            )
            base = (;
                δx = copy(δx),
                δu = copy(δu),
                result,
                iter,
            )
            if opts.keep_first_consistent || opts.objective == :first_consistent
                return BackwardStepRecord(
                    k,
                    :ok,
                    Float64(result.cost),
                    result.E_prev,
                    result.kappa,
                    _record_info(
                        result.L,
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
                        _candidate_diagnostics_empty()...,
                    ),
                )
            end
            break
        end

        opts.verbose && println(
            "[adaptive-box] k=",
            k,
            " iter=",
            iter,
            " inconsistent δx=",
            δx,
            " δu=",
            δu,
            " required_δx=",
            result.required_X_radius,
            " required_δu=",
            result.required_U_radius,
            " L=",
            result.L,
        )

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
                xk,
                xnext,
                uk,
                wk,
                Wbar,
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

            opts.verbose && println(
                "[adaptive-box] k=",
                k,
                " candidate=",
                idx,
                " scale=",
                scale,
                " status=",
                candidate.status,
                " δx=",
                δx_candidate,
                " δu=",
                δu_candidate,
                " required_δx=",
                candidate.required_X_radius,
                " required_δu=",
                candidate.required_U_radius,
                " logvolume=",
                candidate.logvolume,
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
        opts.verbose && println(
            "[adaptive-box] k=",
            k,
            " selected candidate=",
            best_index,
            " scale=",
            best_scale,
            " logvolume=",
            best.result.logvolume,
            " δx=",
            best.δx,
            " δu=",
            best.δu,
        )
        return BackwardStepRecord(
            k,
            :ok,
            Float64(best.result.cost),
            best.result.E_prev,
            best.result.kappa,
            _record_info(
                best.result.L,
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
                diag_tuple...,
            ),
        )
    end

    opts.verbose && println(
        "[adaptive-box] k=",
        k,
        " failed status=",
        last_status,
        " δx=",
        δx,
        " δu=",
        δu,
        " required_δx=",
        last_required_X,
        " required_δu=",
        last_required_U,
        " L=",
        last_L,
    )
    return BackwardStepRecord(
        k,
        :infeasible,
        nothing,
        nothing,
        nothing,
        (;
            L = last_L,
            Xbar_radius = copy(δx),
            Ubar_radius = copy(δu),
            required_X_radius = last_required_X,
            required_U_radius = last_required_U,
            adaptive_box_iters = last_iter,
            adaptive_box_status = last_status,
            selected_logvolume = nothing,
            selected_scale = nothing,
            selected_candidate_index = nothing,
            number_of_candidate_boxes = 0,
            _candidate_diagnostics_empty()...,
        ),
    )
end

function backward_step!(ctx::EllipsoidalBackwardContext, k::Int, E_next)
    ctx.adaptive_boxes.enabled && return _adaptive_backward_step!(ctx, k, E_next)
    return _fixed_backward_step!(ctx, k, E_next)
end

function _collect_kappas(steps::AbstractVector{<:BackwardStepRecord})
    first_idx = findfirst(step -> step.kappa !== nothing, steps)
    first_idx === nothing && return Nothing[]

    κ1 = steps[first_idx].kappa
    kappas = Vector{typeof(κ1)}()

    for step in steps
        κ = step.kappa
        κ === nothing && continue
        push!(kappas, κ)
    end

    return kappas
end

function run_backward_chain!(ctx::EllipsoidalBackwardContext)
    t0 = time()

    nx = length(ctx.xs[end])
    E_next = if ctx.terminal_center !== nothing && ctx.terminal_shape !== nothing
        UT.Ellipsoid(ctx.terminal_shape, ctx.terminal_center)
    else
        PN = Matrix{Float64}(LA.I, nx, nx) * (1.0 / ctx.terminal_radius^2)
        UT.Ellipsoid(PN, collect(ctx.xs[end]))
    end
    

    steps = BackwardStepRecord[]
    ellipsoids = [E_next]

    for k in ctx.K:-1:1
        rec = backward_step!(ctx, k, E_next)
        push!(steps, rec)
       
        if rec.status == :infeasible
            return EllipsoidalCertificationResult(
                false,
                k,
                Float64(time() - t0),
                steps,
                nothing,
                (; ellipsoids, kappas = _collect_kappas(steps)),
            )
        end
        
        #println("[",k,"] is a succes\n")
        E_next = rec.ellipsoid
        
        #println("the volume is ",UT.get_volume(E_next),"\n\n")
        push!(ellipsoids, rec.ellipsoid)
        #println(E_next,"\n\n")
    end

    kappas = _collect_kappas(steps)

    return EllipsoidalCertificationResult(
        true,
        nothing,
        Float64(time() - t0),
        steps,
        kappas,
        (; ellipsoids, kappas),
    )
end
