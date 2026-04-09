export DiscreteSystems

module DiscreteSystems

import Dionysos
const UT = Dionysos.Utils
const ST = Dionysos.System
const PR = Dionysos.Problem

import MathematicalSystems as MS
import JuMP: MOI

# ------------ Control ------------

struct PredicateDomain{F}
    pred::F
end
Base.in(x, X::PredicateDomain) = X.pred(x)

mutable struct SymbolicControlTable
    U::Vector{Vector{Int}}   # U[q] = admissible symbols
end

SymbolicControlTable(nstates::Int) = SymbolicControlTable([Int[] for _ in 1:nstates])

add_control!(C::SymbolicControlTable, q::Int, u::Int) = push!(C.U[q], u)

function set_control!(C::SymbolicControlTable, q::Int, u::Int)
    empty!(C.U[q])
    push!(C.U[q], u)
    return u
end

is_defined(C::SymbolicControlTable, q::Int) = !isempty(C.U[q])
function to_ms_controller(C::SymbolicControlTable)
    qfun = (qs::Int) -> C.U[qs] # set-valued output
    X = PredicateDomain((qs::Int) -> is_defined(C, qs))
    return MS.ConstrainedBlackBoxMap(1, 1, qfun, X)
end

include("optimal_control_problem.jl")
include("safety_problem.jl")
include("cosafe_ltl_problem.jl")

end
