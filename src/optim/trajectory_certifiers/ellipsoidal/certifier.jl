# The AbstractTrajectoryCertifier front. Both directions share every verb —
# fields, resets, result accessors, controller extraction — through
# `AbstractEllipsoidalCertifier`; only the chain runner (`_run_certification`)
# differs per direction.

abstract type AbstractEllipsoidalCertifier <: AbstractTrajectoryCertifier end

"""
    BackwardCertifier(affine_provider, backend, options::ChainOptions)

Ellipsoidal backward trajectory certifier: chains one
`ST.solve_transition_backward` SDP per transition, from the terminal ellipsoid back
to the trajectory start, gating every step for soundness (see [`ChainOptions`](@ref)).
"""
mutable struct BackwardCertifier{AP, Backend, Opts} <: AbstractEllipsoidalCertifier
    problem::Union{Nothing, PR.ProblemType}
    traj::Union{Nothing, ST.Trajectory}

    affine_provider::AP
    backend::Backend
    options::Opts

    result::Union{Nothing, CertificationResult}
    success::Bool
    solve_time_sec::Float64
end

function BackwardCertifier(affine_provider, backend, options::ChainOptions)
    return BackwardCertifier(
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

"""
    ForwardCertifier(affine_provider, backend, options::ForwardOptions)

Ellipsoidal forward trajectory certifier: propagates a certified tube from the
entry ellipsoid along the nominal trajectory (see [`ForwardOptions`](@ref)).
Success requires every step certified, every enabled gate passed, and — when
`check_terminal` — the final tube inside the problem's target set. The result
reports the per-step contraction profile.
"""
mutable struct ForwardCertifier{AP, Backend, Opts} <: AbstractEllipsoidalCertifier
    problem::Union{Nothing, PR.ProblemType}
    traj::Union{Nothing, ST.Trajectory}

    affine_provider::AP
    backend::Backend
    options::Opts

    result::Union{Nothing, CertificationResult}
    success::Bool
    solve_time_sec::Float64
end

function ForwardCertifier(affine_provider, backend, options::ForwardOptions)
    return ForwardCertifier(
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

function _reset!(cert::AbstractEllipsoidalCertifier)
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function set_problem!(cert::AbstractEllipsoidalCertifier, prob::PR.ProblemType)
    cert.problem = prob
    return _reset!(cert)
end

function set_trajectory!(cert::AbstractEllipsoidalCertifier, traj::ST.Trajectory)
    cert.traj = traj
    return _reset!(cert)
end

function certify!(cert::AbstractEllipsoidalCertifier)
    @assert cert.problem !== nothing "Call set_problem!(cert, problem) first."
    @assert cert.traj !== nothing "Call set_trajectory!(cert, trajectory) first."

    ctx = build_context(
        cert.problem,
        cert.traj,
        cert.affine_provider,
        cert.backend,
        cert.options,
    )
    res = _run_certification(cert, ctx)

    cert.result = res
    cert.success = res.success
    cert.solve_time_sec = res.solve_time_sec
    return cert
end

_run_certification(::BackwardCertifier, ctx::ChainContext) = run_chain!(ctx)
_run_certification(::ForwardCertifier, ctx::ChainContext) = _run_forward_chain(ctx)

"""
    get_result(cert) -> Union{Nothing, CertificationResult}

The last [`CertificationResult`](@ref), or `nothing` before `certify!`.
"""
get_result(cert::AbstractEllipsoidalCertifier) = cert.result
get_success(cert::AbstractEllipsoidalCertifier) = cert.success
get_solve_time(cert::AbstractEllipsoidalCertifier) = cert.solve_time_sec

# The certified controller is a real controller (`ST.FunnelController`), usable by
# the closed-loop simulation protocol — not a bare vector of gains.
function get_controller(cert::AbstractEllipsoidalCertifier)
    res = cert.result
    (res === nothing || !res.success) && return nothing
    return ST.FunnelController(
        collect(res.lmi_data.kappas),
        collect(res.lmi_data.ellipsoids),
    )
end
