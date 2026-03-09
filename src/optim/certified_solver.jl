export CertifiedSolveResult, CertifiedPipelineSolver,
       set_problem!, solve!, get_result, get_success, get_solve_time

const SC = SymbolicCertifier

struct CertifiedSolveResult{CAND, CERT}
    candidate::CAND
    certification::CERT
    success::Bool
    solve_time_sec::Float64
end

mutable struct CertifiedPipelineSolver{G, C, P, R}
    generator::G
    certifier::C
    problem::Union{Nothing, P}
    result::Union{Nothing, R}
end

function set_problem!(solver::CertifiedPipelineSolver, prob)
    solver.problem = prob
    solver.result = nothing
    set_problem!(solver.generator, prob)
    SC.set_problem!(solver.certifier, prob)
    return solver
end

function solve!(solver::CertifiedPipelineSolver)
    @assert solver.problem !== nothing "Call set_problem!(solver, problem) first."

    cand = nothing
    certres = nothing
    ok = false

    t = @elapsed begin
        generate!(solver.generator)
        cand = get_trajectory(solver.generator)
        @assert cand !== nothing "Generator returned nothing trajectory."
        SC.set_trajectory!(solver.certifier, cand)
        SC.certify!(solver.certifier)
        certres = SC.get_result(solver.certifier)
        @assert certres !== nothing "Certifier returned nothing result."
        ok = get_success(solver.generator) && SC.get_success(solver.certifier)
    end

    solver.result = CertifiedSolveResult(cand, certres, ok, t)
    return solver
end

get_result(solver::CertifiedPipelineSolver) = solver.result
get_success(solver::CertifiedPipelineSolver) = solver.result === nothing ? false : solver.result.success
get_solve_time(solver::CertifiedPipelineSolver) = solver.result === nothing ? NaN : solver.result.solve_time_sec