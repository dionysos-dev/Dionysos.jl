
export AbstractTrajectoryCertifier,
    set_problem!, set_trajectory!, certify!, get_controller, get_success, get_solve_time

abstract type AbstractTrajectoryCertifier end

function set_problem!(cert::AbstractTrajectoryCertifier, concrete_problem)
    return error("set_problem! not implemented for $(typeof(cert))")
end

function set_trajectory!(cert::AbstractTrajectoryCertifier, traj)
    return error("set_trajectory! not implemented for $(typeof(cert))")
end

function certify!(cert::AbstractTrajectoryCertifier)
    return error("certify! not implemented for $(typeof(cert))")
end

get_controller(cert::AbstractTrajectoryCertifier) = nothing
get_success(cert::AbstractTrajectoryCertifier) = false
get_solve_time(cert::AbstractTrajectoryCertifier) = NaN

include("uniform_grid_trajectory_certifier.jl")
include("ellipsoidal/ellipsoidal.jl")
