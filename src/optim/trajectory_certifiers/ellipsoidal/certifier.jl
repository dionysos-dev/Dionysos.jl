# The AbstractTrajectoryCertifier front: holds problem + trajectory + provider +
# backend + options; `certify!` builds the context and runs the chain.

"""
    BackwardCertifier(affine_provider, backend, options::ChainOptions)

Ellipsoidal backward trajectory certifier: chains one
`ST.solve_transition_backward` SDP per transition, from the terminal ellipsoid back
to the trajectory start, gating every step for soundness (see [`ChainOptions`](@ref)).
"""
mutable struct BackwardCertifier{AP, Backend, Opts} <: AbstractTrajectoryCertifier
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

function set_problem!(cert::BackwardCertifier, prob::PR.ProblemType)
    cert.problem = prob
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function set_trajectory!(cert::BackwardCertifier, traj::ST.Trajectory)
    cert.traj = traj
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function certify!(cert::BackwardCertifier)
    @assert cert.problem !== nothing "Call set_problem!(cert, problem) first."
    @assert cert.traj !== nothing "Call set_trajectory!(cert, trajectory) first."

    t0 = time()

    ctx = build_context(
        cert.problem,
        cert.traj,
        cert.affine_provider,
        cert.backend,
        cert.options,
    )

    res = run_chain!(ctx)

    cert.result = res
    cert.success = res.success
    cert.solve_time_sec = time() - t0

    return cert
end

"""
    get_result(cert::BackwardCertifier) -> Union{Nothing, CertificationResult}

The last [`CertificationResult`](@ref), or `nothing` before `certify!`.
"""
get_result(cert::BackwardCertifier) = cert.result
get_success(cert::BackwardCertifier) = cert.success
get_solve_time(cert::BackwardCertifier) = cert.solve_time_sec

get_controller(cert::BackwardCertifier) =
    cert.result === nothing ? nothing : cert.result.controller
