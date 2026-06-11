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

struct EllipsoidalBackwardOptions
    maxδx::Float64
    maxδu::Float64
    λ::Float64
    terminal_radius::Float64
    state_scaling::Union{Nothing, Matrix{Float64}}
    terminal_center::Union{Nothing, Vector{Float64}}
    terminal_shape::Union{Nothing, Matrix{Float64}}
    adaptive_boxes::AdaptiveLinearizationBoxOptions
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
    options::EllipsoidalBackwardOptions
end

function _identity_transition_cost(nx::Int, nu::Int)
    return [
        Matrix{Float64}(LA.I, nx, nx) zeros(nx, nu) zeros(nx, 1)
        zeros(nu, nx) Matrix{Float64}(LA.I, nu, nu) zeros(nu, 1)
        zeros(1, nx + nu) 1.0
    ]
end

function build_symbolic_context(
    problem::DI.Problem.ProblemType,
    candidate::DI.Optim.CandidateTrajectory,
    config::EllipsoidalBackwardConfig,
    symbolic_builder,
)
    xs = collect(ST.enum_elems(candidate.x_traj))
    us = collect(ST.enum_elems(candidate.u_traj))
    K = length(us)
    @assert length(xs) == K + 1

    sym = symbolic_builder(problem, candidate, config)

    nx = length(xs[1])
    nu = length(us[1])
    S = _identity_transition_cost(nx, nu)

    opts = EllipsoidalBackwardOptions(config.options, sym, nx, nu)

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
        opts,
    )
end

_floatvec(x) = Float64.(collect(x))

function _interval_box_radii(Δ)
    return [max(abs(Float64(IA.inf(I))), abs(Float64(IA.sup(I)))) for I in collect(Δ)]
end

function AdaptiveLinearizationBoxOptions(opts, sym, nx::Int, nu::Int)
    default_ΔX = _interval_box_radii(sym.ΔX)
    default_ΔU = _interval_box_radii(sym.ΔU)

    ΔX_min = _floatvec(get(opts, :ΔX_min, zeros(nx)))
    ΔU_min = _floatvec(get(opts, :ΔU_min, zeros(nu)))

    return AdaptiveLinearizationBoxOptions(
        Bool(get(opts, :adaptive_linearization_boxes, false)),
        max.(_floatvec(get(opts, :ΔX_initial, default_ΔX)), ΔX_min),
        ΔX_min,
        _floatvec(get(opts, :ΔX_max, fill(Inf, nx))),
        max.(_floatvec(get(opts, :ΔU_initial, default_ΔU)), ΔU_min),
        ΔU_min,
        _floatvec(get(opts, :ΔU_max, fill(Inf, nu))),
        Float64(get(opts, :adaptive_box_growth, 2.0)),
        Float64(get(opts, :adaptive_box_safety, 1.05)),
        Int(get(opts, :adaptive_box_max_iters, 1)),
        Float64(get(opts, :adaptive_box_atol, 1.0e-8)),
        Bool(get(opts, :adaptive_box_verbose, false)),
        _floatvec(
            get(opts, :adaptive_box_search_scales, [0.7, 0.85, 1.0, 1.15, 1.35, 1.6, 2.0]),
        ),
        Symbol(get(opts, :adaptive_box_objective, :max_volume)),
        Bool(get(opts, :adaptive_box_keep_first_consistent, false)),
    )
end

function EllipsoidalBackwardOptions(opts, sym, nx::Int, nu::Int)
    terminal_center, terminal_shape = _terminal_ellipsoid_data(
        get(opts, :terminal_center, nothing),
        get(opts, :terminal_shape, nothing),
    )

    return EllipsoidalBackwardOptions(
        Float64(opts.maxδx),
        Float64(opts.maxδu),
        Float64(opts.λ),
        Float64(opts.rayon_terminal),
        _state_scaling(get(opts, :state_scaling, nothing)),
        terminal_center,
        terminal_shape,
        AdaptiveLinearizationBoxOptions(opts, sym, nx, nu),
    )
end

EllipsoidalBackwardOptions(opts::EllipsoidalBackwardOptions, sym, nx::Int, nu::Int) = opts

function _state_scaling(raw)
    raw === nothing && return nothing

    if raw isa AbstractMatrix
        return Matrix{Float64}(raw)
    end

    return Matrix{Float64}(LA.Diagonal(_floatvec(raw)))
end

function _terminal_ellipsoid_data(raw_center, raw_shape)
    raw_center === nothing && return nothing, nothing
    return _floatvec(raw_center), Matrix(LA.Symmetric(Matrix{Float64}(raw_shape)))
end

function _scale_target_ellipsoid(E, xnext, scaling)
    scaling === nothing && return E
    # Convention used here: E(c, P) = {x : (x-c)' P (x-c) <= 1}.
    # With x = xnext + T*z, the normalized ellipsoid center is
    # cz = T^{-1} * (c - xnext), not necessarily zero.
    T = Matrix{Float64}(scaling)
    Tinv = inv(T)
    Pz = Matrix(T' * Matrix(E.P) * T)
    cz = collect(Tinv * (collect(E.c) - collect(xnext)))
    return UT.Ellipsoid(Pz, cz)
end

function _unscale_source_ellipsoid(Ez, center, scaling)
    scaling === nothing && return Ez
    # Convention used here: E(c, P) = {x : (x-c)' P (x-c) <= 1}.
    # Since z = T^{-1}(x-c), the physical ellipsoid has Px = T^{-T} * Pz * T^{-1}.
    Tinv = inv(Matrix{Float64}(scaling))
    Px = Matrix(Tinv' * Matrix(Ez.P) * Tinv)
    return UT.Ellipsoid(Px, center)
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

function _scaled_lipschitz(L, nx::Int, scaling)
    scaling === nothing && return L
    Ls = collect(Float64, L)
    # The first nx entries are interpreted as output-state remainder radii.
    # For z_next = T^{-1}(x_next - xnext_nominal), bound the transformed
    # interval remainder by the axis-aligned enclosure abs(T^{-1}) * Lx.
    Tinv = inv(Matrix{Float64}(scaling))
    Ls[1:nx] .= abs.(Tinv) * Ls[1:nx]
    return Ls
end

function _scaled_transition_cost(S, xk, nu::Int, scaling)
    scaling === nothing && return S

    # In transition_backward, S is used as a linear factor S * [x; u; 1].
    # With x = xk + Tstate*z, use [x; u; 1] = T * [z; u; 1], hence Sz = S*T.
    nx = length(xk)
    T = zeros(nx + nu + 1, nx + nu + 1)
    T[1:nx, 1:nx] .= Matrix{Float64}(scaling)
    T[1:nx, end] .= xk
    T[(nx + 1):(nx + nu), (nx + 1):(nx + nu)] .= Matrix{Float64}(LA.I, nu, nu)
    T[end, end] = 1.0
    return S * T
end

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
    Kx = Kz * inv(Matrix{Float64}(state_scaling))
    cont = MS.AffineMap(Kx, ell - Kx * xk)
    return E_prev, cont, cost
end

function _centered_interval_box(center::AbstractVector, radius::AbstractVector)
    c = Float64.(collect(center))
    r = Float64.(collect(radius))
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

function _controller_matrices(kappa::MS.AffineMap, nx::Int)
    return Matrix{Float64}(kappa.A), vec(Float64.(kappa.c))
end

function _controller_matrices(kappa::AbstractMatrix, nx::Int)
    K = Matrix{Float64}(kappa[:, 1:nx])
    b = vec(Float64.(kappa[:, nx + 1]))
    return K, b
end

function _controller_matrices(kappa, nx::Int)
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

_grow_infeasible_box_radii(δ::AbstractVector, δmax::AbstractVector, growth::Float64) =
    min.(growth .* δ, δmax)

_grow_to_required_box_radii(
    required::AbstractVector,
    δmin::AbstractVector,
    δmax::AbstractVector,
    safety::Float64,
) = min.(max.(safety .* required, δmin), δmax)

function _solve_transition(
    ctx::EllipsoidalBackwardContext,
    affineSys,
    E_next,
    xk,
    xnext,
    uk,
    L,
)
    return _transition_backward(
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
        maxδx = ctx.options.maxδx,
        maxδu = ctx.options.maxδu,
        λ = ctx.options.λ,
        state_scaling = ctx.options.state_scaling,
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

    affineSys, L = ST.buildAffineApproximation(
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

    E_prev, kappa, cost = _solve_transition(ctx, affineSys, E_next, xk, xnext, uk, L)
    Xbar_radius = _interval_box_radii(ctx.symbolic.ΔX)
    Ubar_radius = _interval_box_radii(ctx.symbolic.ΔU)

    if E_prev === nothing || kappa === nothing
        return BackwardStepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            _step_summary(
                L,
                Xbar_radius,
                Ubar_radius,
                nothing,
                nothing,
                1,
                :fixed_infeasible;
                _candidate_diagnostics_empty()...,
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
        _step_summary(
            L,
            Xbar_radius,
            Ubar_radius,
            required_X_radius,
            required_U_radius,
            1,
            :fixed;
            _candidate_diagnostics_empty()...,
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

function _step_summary(
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
    affineSys, L = ST.buildAffineApproximation(
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
    E_prev, kappa, cost = _solve_transition(ctx, affineSys, E_next, xk, xnext, uk, L)

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

    logvolume = _ellipsoid_logvolume(E_prev)
    if !isfinite(logvolume)
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
    opts = ctx.options.adaptive_boxes
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

        if result.status == :ok
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
            base = (; δx = copy(δx), δu = copy(δu), result, iter)
            if opts.keep_first_consistent || opts.objective == :first_consistent
                return BackwardStepRecord(
                    k,
                    :ok,
                    Float64(result.cost),
                    result.E_prev,
                    result.kappa,
                    _step_summary(
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
            _step_summary(
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
        _step_summary(
            last_L,
            δx,
            δu,
            last_required_X,
            last_required_U,
            last_iter,
            last_status;
            _candidate_diagnostics_empty()...,
        ),
    )
end

function backward_step!(ctx::EllipsoidalBackwardContext, k::Int, E_next)
    ctx.options.adaptive_boxes.enabled && return _adaptive_backward_step!(ctx, k, E_next)
    return _fixed_backward_step!(ctx, k, E_next)
end

function _collect_kappas(steps::AbstractVector{<:BackwardStepRecord})
    valid_steps = filter(step -> step.kappa !== nothing, steps)
    isempty(valid_steps) && return Nothing[]
    return getproperty.(valid_steps, :kappa)
end

function run_backward_chain!(ctx::EllipsoidalBackwardContext)
    t0 = time()

    nx = length(ctx.xs[end])
    opts = ctx.options
    E_next = if opts.terminal_center !== nothing && opts.terminal_shape !== nothing
        UT.Ellipsoid(opts.terminal_shape, opts.terminal_center)
    else
        PN = Matrix{Float64}(LA.I, nx, nx) * (1.0 / opts.terminal_radius^2)
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

        E_next = rec.ellipsoid
        push!(ellipsoids, rec.ellipsoid)
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
