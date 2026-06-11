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

mutable struct EllipsoidalBackwardCertifier{CFG, SB} <: AbstractSymbolicCertifier
    problem::Union{Nothing, Dionysos.Problem.ProblemType}
    candidate::Union{Nothing, Dionysos.Optim.CandidateTrajectory}
    config::CFG
    result::Union{Nothing, EllipsoidalCertificationResult}
    success::Bool
    solve_time_sec::Float64
    symbolic_builder::SB
end

function EllipsoidalBackwardCertifier(config::EllipsoidalBackwardConfig, symbolic_builder)
    return EllipsoidalBackwardCertifier(
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

function set_trajectory!(
    cert::EllipsoidalBackwardCertifier,
    cand::Dionysos.Optim.CandidateTrajectory,
)
    cert.candidate = cand
    cert.result = nothing
    cert.success = false
    cert.solve_time_sec = 0.0
    return cert
end

function certify!(cert::EllipsoidalBackwardCertifier)
    @assert cert.problem !== nothing "Call set_problem!(cert, problem) first."
    @assert cert.candidate !== nothing "Call set_trajectory!(cert, candidate) first."

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
