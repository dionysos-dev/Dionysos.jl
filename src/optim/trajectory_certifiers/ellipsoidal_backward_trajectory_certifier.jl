module EllipsoidalBackwardTrajectoryCertifier

import Dionysos
import LinearAlgebra as LA
import MathematicalSystems as MS

const DI = Dionysos
const PR = DI.Problem
const UT = DI.Utils
const ST = DI.System

import ..AbstractTrajectoryCertifier
import ..set_problem!
import ..set_trajectory!
import ..certify!
import ..get_controller
import ..get_success
import ..get_solve_time

export TrajectoryCertifier,
    AdaptiveLinearizationBoxOptions,
    EllipsoidalBackwardOptions,
    BackwardStepRecord,
    EllipsoidalCertificationResult,
    get_result,
    build_symbolic_context,
    backward_step!,
    run_backward_chain!

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

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

Base.@kwdef mutable struct EllipsoidalBackwardOptions
    maxδx::Float64 = 0.2
    maxδu::Float64 = 0.5
    λ::Float64 = 1.0
    terminal_shape::Union{Nothing, Matrix{Float64}} = nothing
    transition_cost::Union{Nothing, UT.ScalarControlFunction, Matrix{Float64}} = nothing
    state_scaling::Union{Nothing, Matrix{Float64}} = nothing

    # Used when adaptive_boxes === nothing or adaptive_boxes.enabled == false
    linearization_δx::Vector{Float64} = Float64[]
    linearization_δu::Vector{Float64} = Float64[]

    adaptive_boxes::Union{Nothing, AdaptiveLinearizationBoxOptions} = nothing
    use_log_det::Bool = true
end

# ------------------------------------------------------------
# Result types
# ------------------------------------------------------------

struct BackwardStepRecord{TE, TK, TS}
    k::Int
    status::Symbol
    cost::Union{Nothing, Float64}
    ellipsoid::TE
    kappa::TK
    summary::TS
end

function BackwardStepRecord(
    k::Integer,
    status::Symbol,
    cost,
    ellipsoid,
    kappa,
    summary::TS,
) where {TS}
    return BackwardStepRecord(
        Int(k),
        status,
        cost === nothing ? nothing : Float64(cost),
        ellipsoid,
        kappa,
        summary,
    )
end

struct EllipsoidalCertificationResult{S, CTRL, LMI}
    success::Bool
    failed_k::Union{Nothing, Int}
    solve_time_sec::Float64
    steps::Vector{S}
    controller::CTRL
    lmi_data::LMI
end

function EllipsoidalCertificationResult(
    success::Bool,
    failed_k::Union{Nothing, Integer},
    solve_time_sec::Real,
    steps::Vector{S},
    controller,
    lmi_data,
) where {S}
    return EllipsoidalCertificationResult(
        success,
        failed_k === nothing ? nothing : Int(failed_k),
        Float64(solve_time_sec),
        steps,
        controller,
        lmi_data,
    )
end

# ------------------------------------------------------------
# Certifier
# ------------------------------------------------------------

mutable struct TrajectoryCertifier{AP, Backend, Opts} <: AbstractTrajectoryCertifier
    problem::Union{Nothing, PR.ProblemType}
    traj::Union{Nothing, ST.ClosedLoopTrajectory}

    affine_provider::AP
    backend::Backend
    options::Opts

    result::Union{Nothing, EllipsoidalCertificationResult}
    success::Bool
    solve_time_sec::Float64
end

function TrajectoryCertifier(affine_provider, backend, options::EllipsoidalBackwardOptions)
    return TrajectoryCertifier(
        nothing,
        nothing,
        affine_provider,
        backend,
        options,
        nothing,
        false,
        0.0,
    )
end

function set_problem!(cert::TrajectoryCertifier, prob::PR.ProblemType)
    cert.problem = prob
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function set_trajectory!(cert::TrajectoryCertifier, traj::ST.ClosedLoopTrajectory)
    cert.traj = traj
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function certify!(cert::TrajectoryCertifier)
    @assert cert.problem !== nothing "Call set_problem!(cert, problem) first."
    @assert cert.traj !== nothing "Call set_trajectory!(cert, trajectory) first."

    t0 = time()

    ctx = build_symbolic_context(
        cert.problem,
        cert.traj,
        cert.affine_provider,
        cert.backend,
        cert.options,
    )

    res = run_backward_chain!(ctx)

    cert.result = res
    cert.success = res.success
    cert.solve_time_sec = time() - t0

    return cert
end

get_result(cert::TrajectoryCertifier) = cert.result
get_success(cert::TrajectoryCertifier) = cert.success
get_solve_time(cert::TrajectoryCertifier) = cert.solve_time_sec

get_controller(cert::TrajectoryCertifier) =
    cert.result === nothing ? nothing : cert.result.controller

# ------------------------------------------------------------
# Context
# ------------------------------------------------------------

struct EllipsoidalBackwardContext{P, TR, AP, TX, TU, TS, TB, OPTS}
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
    if cost === nothing
        return _identity_transition_cost(nx, nu)
    elseif cost isa Matrix
        S = Matrix{Float64}(cost)
    elseif cost isa UT.QuadraticStateControlFunction
        S = Matrix{Float64}(UT.get_full_psd_matrix(cost))
    else
        error("Unsupported transition cost type: $(typeof(cost))")
    end

    @assert size(S) == (nx + nu + 1, nx + nu + 1)
    return S
end

function build_symbolic_context(
    problem::PR.ProblemType,
    traj::ST.ClosedLoopTrajectory,
    affine_provider,
    backend,
    options::EllipsoidalBackwardOptions,
)
    xs = collect(ST.enum_elems(traj.x))
    us = collect(ST.enum_elems(traj.u))

    K = length(us)
    @assert length(xs) == K + 1 "Expected length(xs) == length(us) + 1."

    nx = length(xs[1])
    nu = length(us[1])
    S = _transition_cost_matrix(options.transition_cost, nx, nu)

    return EllipsoidalBackwardContext(
        problem,
        traj,
        affine_provider,
        xs,
        us,
        K,
        S,
        backend,
        options,
    )
end

# ------------------------------------------------------------
# Scaling helpers
# ------------------------------------------------------------

function _scale_target_ellipsoid(E, xnext, scaling)
    scaling === nothing && return E

    T = Matrix{Float64}(scaling)
    Tinv = inv(T)
    Pz = Matrix(T' * Matrix(E.P) * T)
    cz = collect(Tinv * (collect(E.c) - collect(xnext)))

    return UT.Ellipsoid(Pz, cz)
end

function _unscale_source_ellipsoid(Ez, center, scaling)
    scaling === nothing && return Ez

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
    Tinv = inv(Matrix{Float64}(scaling))
    Ls[1:nx] .= abs.(Tinv) * Ls[1:nx]

    return Ls
end

function _scaled_transition_cost(S, xk, nu::Int, scaling)
    scaling === nothing && return S

    nx = length(xk)
    T = zeros(nx + nu + 1, nx + nu + 1)

    T[1:nx, 1:nx] .= Matrix{Float64}(scaling)
    T[1:nx, end] .= xk
    T[(nx + 1):(nx + nu), (nx + 1):(nx + nu)] .= Matrix{Float64}(LA.I, nu, nu)
    T[end, end] = 1.0

    return S * T
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
    use_log_det = true,
)
    if state_scaling === nothing
        return UT.transition_backward(
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
            use_log_det = use_log_det,
        )
    end

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

function _solve_transition(ctx::EllipsoidalBackwardContext, approx, E_next, xk, xnext, uk)
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
        use_log_det = ctx.options.use_log_det,
    )
end

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

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

function _box_contains_required_radii(δx, δu, required_δx, required_δu; atol::Float64)
    return all(required_δx .<= δx .+ atol) && all(required_δu .<= δu .+ atol)
end

_clamp_box_radii(δ, δmin, δmax) = min.(max.(δ, δmin), δmax)
_grow_infeasible_box_radii(δ, δmax, growth) = min.(growth .* δ, δmax)
_grow_to_required_box_radii(required, δmin, δmax, safety) =
    min.(max.(safety .* required, δmin), δmax)

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
    provider_summary = nothing,
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
        provider_summary,
    )
end

# ------------------------------------------------------------
# Backward steps
# ------------------------------------------------------------

function _fixed_backward_step!(ctx::EllipsoidalBackwardContext, k::Int, E_next)
    opts = ctx.options

    @assert !isempty(opts.linearization_δx) "Set options.linearization_δx for fixed mode."
    @assert !isempty(opts.linearization_δu) "Set options.linearization_δu for fixed mode."

    xk = collect(ctx.xs[k])
    xnext = collect(E_next.c)
    uk = collect(ctx.us[k])

    δx = collect(Float64, opts.linearization_δx)
    δu = collect(Float64, opts.linearization_δu)

    skip_already_in = false #make this a param solver
    if skip_already_in && xk ∈ E_next
        return BackwardStepRecord(
            k,
            :already_inside,
            0.0,
            E_next,
            nothing,
            _step_summary(
                nothing,
                δx,
                δu,
                nothing,
                nothing,
                0,
                :already_inside_next_ellipsoid;
                provider_summary = nothing,
                _candidate_diagnostics_empty()...,
            ),
        )
    end

    u_lin = uk # zeros(length(uk))
    approx = ST.build_affine_approximation(ctx.affine_provider, k, xk, xnext, u_lin, δx, δu)

    E_prev, kappa, cost = _solve_transition(ctx, approx, E_next, xk, xnext, uk)

    if E_prev === nothing || kappa === nothing
        return BackwardStepRecord(
            k,
            :infeasible,
            nothing,
            nothing,
            nothing,
            _step_summary(
                approx.lipschitz,
                δx,
                δu,
                nothing,
                nothing,
                1,
                :fixed_infeasible;
                provider_summary = approx.summary,
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
            approx.lipschitz,
            δx,
            δu,
            required_X_radius,
            required_U_radius,
            1,
            :fixed;
            provider_summary = approx.summary,
            _candidate_diagnostics_empty()...,
        ),
    )
end

function _evaluate_adaptive_box_candidate(ctx, E_next, k, xk, xnext, uk, δx, δu; atol)
    u_lin = uk # zeros(length(uk))
    approx = ST.build_affine_approximation(ctx.affine_provider, k, xk, xnext, u_lin, δx, δu)

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
    @assert opts !== nothing "adaptive_boxes cannot be nothing in adaptive mode."

    xk = collect(ctx.xs[k])
    xnext = collect(E_next.c)
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
                return BackwardStepRecord(
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

        return BackwardStepRecord(
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

    return BackwardStepRecord(
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

function backward_step!(ctx::EllipsoidalBackwardContext, k::Int, E_next)
    opts = ctx.options.adaptive_boxes
    if opts !== nothing && opts.enabled
        return _adaptive_backward_step!(ctx, k, E_next)
    end

    return _fixed_backward_step!(ctx, k, E_next)
end

# ------------------------------------------------------------
# Backward chain
# ------------------------------------------------------------

function _collect_kappas(steps::AbstractVector{<:BackwardStepRecord})
    valid_steps = filter(step -> step.kappa !== nothing, steps)
    isempty(valid_steps) && return Nothing[]
    return getproperty.(valid_steps, :kappa)
end

function run_backward_chain!(ctx::EllipsoidalBackwardContext)
    t0 = time()

    nx = length(ctx.xs[end])
    opts = ctx.options

    @assert opts.terminal_shape !== nothing "Set options.terminal_shape."

    terminal_shape = Matrix{Float64}(opts.terminal_shape)

    @assert size(terminal_shape) == (nx, nx) "
    terminal_shape must have size ($nx, $nx).
    "

    E_next = UT.Ellipsoid(terminal_shape, collect(ctx.xs[end]))

    steps = BackwardStepRecord[]
    ellipsoids = [E_next]

    k = ctx.K

    while k >= 1
        rec = backward_step!(ctx, k, E_next)

        skip_point = false
        if skip_point && rec.status == :infeasible && k > 1
            # Skip y_k and try y_{k-1} -> current E_next,
            # where E_next is still centered around y_{k+1}.
            rec_skip = backward_step!(ctx, k - 1, E_next)

            if rec_skip.status != :infeasible
                rec = BackwardStepRecord(
                    rec_skip.k,
                    :ok_skip,
                    rec_skip.cost,
                    rec_skip.ellipsoid,
                    rec_skip.kappa,
                    merge(
                        rec_skip.summary,
                        (;
                            skipped = true,
                            skipped_index = k,
                            skipped_state = collect(ctx.xs[k]),
                            reconnected_from = k - 1,
                            reconnected_to = min(k + 1, length(ctx.xs)),
                        ),
                    ),
                )

                push!(steps, rec)

                E_next = rec.ellipsoid
                push!(ellipsoids, rec.ellipsoid)

                k -= 2
                continue
            end
        end

        push!(steps, rec)

        if rec.status == :infeasible
            steps_forward = reverse(steps)
            ellipsoids_forward = reverse(ellipsoids)
            kappas = _collect_kappas(steps_forward)

            return EllipsoidalCertificationResult(
                false,
                k,
                Float64(time() - t0),
                steps_forward,
                nothing,
                (; ellipsoids = ellipsoids_forward, kappas),
            )
        end

        E_next = rec.ellipsoid
        push!(ellipsoids, rec.ellipsoid)

        k -= 1
    end

    steps_forward = reverse(steps)
    ellipsoids_forward = reverse(ellipsoids)
    kappas = _collect_kappas(steps_forward)

    return EllipsoidalCertificationResult(
        true,
        nothing,
        Float64(time() - t0),
        steps_forward,
        kappas,
        (; ellipsoids = ellipsoids_forward, kappas),
    )
end

end # module
