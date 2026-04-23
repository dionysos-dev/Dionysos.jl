module DiscreteSystems

import Dionysos
const ST = Dionysos.System
const PR = Dionysos.Problem

import JuMP: MOI

# ------ Internal data structure for Discrete controllers ------
mutable struct DiscreteControlTable
    U::Vector{Vector{Int}}
end

DiscreteControlTable(nstates::Int) = DiscreteControlTable([Int[] for _ in 1:nstates])

(C::DiscreteControlTable)(q::Int) = C.U[q]

add_control!(C::DiscreteControlTable, q::Int, u::Int) = push!(C.U[q], u)

function set_control!(C::DiscreteControlTable, q::Int, u::Int)
    empty!(C.U[q])
    push!(C.U[q], u)
    return u
end

include("optimal_control_problem.jl")
include("safety_problem.jl")
include("cosafe_ltl_problem.jl")

end
