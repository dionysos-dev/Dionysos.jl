export CertifiedSolveResult, CertifiedPipelineSolver,
       set_problem!, solve!, get_result, get_success, get_solve_time, get_certification_candidate

import Dionysos

const SC = SymbolicCertifier

struct CertifiedSolveResult{CAND, CERTCAND, CERT}
    candidate::CAND
    certification_candidate::CERTCAND
    certification::CERT
    success::Bool
    solve_time_sec::Float64
end

function CertifiedSolveResult(
    candidate,
    certification_candidate,
    certification,
    success::Bool,
    solve_time_sec::Real,
)
    return CertifiedSolveResult(
        candidate,
        certification_candidate,
        certification,
        success,
        Float64(solve_time_sec),
    )
end

mutable struct CertifiedPipelineSolver{G, C, P, R}
    generator::G
    certifier::C
    problem::Union{Nothing, P}
    result::Union{Nothing, R}
end

function CertifiedPipelineSolver(generator, certifier)
    return CertifiedPipelineSolver{
        typeof(generator),
        typeof(certifier),
        Dionysos.Problem.ProblemType,
        CertifiedSolveResult,
    }(
        generator,
        certifier,
        nothing,
        nothing,
    )
end

function set_problem!(solver::CertifiedPipelineSolver, prob)
    solver.problem = prob
    solver.result = nothing
    set_problem!(solver.generator, prob)
    SC.set_problem!(solver.certifier, prob)
    return solver
end

function solve!(solver::CertifiedPipelineSolver; prepare_for_certification = identity)
    @assert solver.problem !== nothing "Call set_problem!(solver, problem) first."

    cand = nothing
    cert_cand = nothing
    certres = nothing
    ok = false

    t = @elapsed begin
        generate!(solver.generator)
        println("we have lunch generate")
        cand = get_trajectory(solver.generator)
        println("we have got our candidate traj")
        @assert cand !== nothing "Generator returned nothing trajectory."

        cert_cand = prepare_for_certification(cand)
        @assert cert_cand !== nothing "prepare_for_certification returned nothing."

        SC.set_trajectory!(solver.certifier, cert_cand)
        SC.certify!(solver.certifier)
        println("we've certify our traj")
        certres = SC.get_result(solver.certifier)
        @assert certres !== nothing "Certifier returned nothing result."

        ok = get_success(solver.generator) && SC.get_success(solver.certifier)
    end

    solver.result = CertifiedSolveResult(cand, cert_cand, certres, ok, t)
    return solver
end

get_result(solver::CertifiedPipelineSolver) = solver.result
get_success(solver::CertifiedPipelineSolver) = solver.result === nothing ? false : solver.result.success
get_solve_time(solver::CertifiedPipelineSolver) = solver.result === nothing ? NaN : solver.result.solve_time_sec
get_certification_candidate(solver::CertifiedPipelineSolver) =
    solver.result === nothing ? nothing : solver.result.certification_candidate