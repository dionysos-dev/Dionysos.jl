import ..SymbolicCertifier:
    AbstractSymbolicCertifier,
    set_problem!,
    set_trajectory!,
    certify!,
    get_success,
    get_solve_time,
    build_symbolic_context,
    run_backward_chain!
import Dionysos

export EllipsoidalBackwardCertifier, get_result

mutable struct EllipsoidalBackwardCertifier{P, C, CFG, R, SB} <: AbstractSymbolicCertifier
    problem::Union{Nothing, P}
    candidate::Union{Nothing, C}
    config::CFG
    result::Union{Nothing, R}
    success::Bool
    solve_time_sec::Float64
    symbolic_builder::SB
end

function EllipsoidalBackwardCertifier(config, symbolic_builder)
    return EllipsoidalBackwardCertifier{
        Dionysos.Problem.ProblemType,
        Dionysos.Optim.CandidateTrajectory,
        typeof(config),
        EllipsoidalCertificationResult,
        typeof(symbolic_builder),
    }(
        nothing,
        nothing,
        config,
        nothing,
        false,
        0.0,
        symbolic_builder,
    )
end

function set_problem!(
    cert::EllipsoidalBackwardCertifier,
    prob::Dionysos.Problem.ProblemType,
)
    cert.problem = prob
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function set_trajectory!(cert::EllipsoidalBackwardCertifier, cand)
    cert.candidate = cand
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function certify!(cert::EllipsoidalBackwardCertifier)
    t0 = time()

    ctx = build_symbolic_context(
        cert.problem,
        cert.candidate,
        cert.config,
        cert.symbolic_builder,
    )

    res = run_backward_chain!(ctx)

    cert.result = res
    cert.success = res.success
    cert.solve_time_sec = time() - t0
    return cert
end

get_result(cert::EllipsoidalBackwardCertifier) = cert.result
get_success(cert::EllipsoidalBackwardCertifier) = cert.success
get_solve_time(cert::EllipsoidalBackwardCertifier) = cert.solve_time_sec
