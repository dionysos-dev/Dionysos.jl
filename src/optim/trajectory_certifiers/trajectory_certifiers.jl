
export AbstractTrajectoryCertifier,
    set_problem!,
    set_trajectory!,
    certify!,
    get_controller,
    get_success,
    get_solve_time,
    get_result

abstract type AbstractTrajectoryCertifier end

# Interface stubs error loudly: a missing method must be a MethodError-grade bug,
# not a silent `nothing`/`false` wrong answer (conventions §1).
function set_problem!(cert::AbstractTrajectoryCertifier, concrete_problem)
    return error("set_problem! not implemented for $(typeof(cert))")
end

function set_trajectory!(cert::AbstractTrajectoryCertifier, traj)
    return error("set_trajectory! not implemented for $(typeof(cert))")
end

function certify!(cert::AbstractTrajectoryCertifier)
    return error("certify! not implemented for $(typeof(cert))")
end

function get_controller(cert::AbstractTrajectoryCertifier)
    return error("get_controller not implemented for $(typeof(cert))")
end

function get_success(cert::AbstractTrajectoryCertifier)
    return error("get_success not implemented for $(typeof(cert))")
end

function get_solve_time(cert::AbstractTrajectoryCertifier)
    return error("get_solve_time not implemented for $(typeof(cert))")
end

"""
    get_result(cert::AbstractTrajectoryCertifier)

The certifier's last result object (certifier-specific type), or `nothing`
before `certify!`. Part of the certifier interface — the bidirectional handoff
and the re-planning loop consume it generically.
"""
function get_result(cert::AbstractTrajectoryCertifier)
    return error("get_result not implemented for $(typeof(cert))")
end

include("uniform_grid_trajectory_certifier.jl")
include("ellipsoidal/ellipsoidal.jl")
