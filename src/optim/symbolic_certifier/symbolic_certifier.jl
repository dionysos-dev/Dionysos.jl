
module SymbolicCertifier

export AbstractSymbolicCertifier,
    set_problem!, set_trajectory!, certify!, get_controller, get_success, get_solve_time
export EllipsoidalBackwardConfig, BackwardStepRecord, EllipsoidalCertificationResult
export EllipsoidalBackwardCertifier

import Dionysos
const PR = Dionysos.Problem

abstract type AbstractSymbolicCertifier end

function set_problem!(cert::AbstractSymbolicCertifier, concrete_problem::PR.ProblemType)
    return error("set_problem! not implemented for $(typeof(cert))")
end

function set_trajectory!(cert::AbstractSymbolicCertifier, traj)
    return error("set_trajectory! not implemented for $(typeof(cert))")
end

function certify!(cert::AbstractSymbolicCertifier)
    return error("certify! not implemented for $(typeof(cert))")
end

get_controller(cert::AbstractSymbolicCertifier) = nothing
get_success(cert::AbstractSymbolicCertifier) = false
get_solve_time(cert::AbstractSymbolicCertifier) = NaN

include("ellipsoidal_backward_types.jl")
include("ellipsoidal_backward_core.jl")
include("ellipsoidal_backward_certifier.jl")
include("uniform_grid_local_tube.jl")

end # SymbolicCertifier module
