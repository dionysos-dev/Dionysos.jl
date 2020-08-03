module Abstraction

using LinearAlgebra
using ProgressMeter

@enum INCL_MODE INNER OUTER

abstract type Automaton end

mutable struct AutomatonList <: Automaton
    nstates::Int
    nsymbols::Int
    transitions::Vector{Tuple{Int, Int, Int}}
    issorted::Bool
end

include("rectangle.jl")
include("gridspace.jl")
include("controlsystem.jl")
include("symbolicmodel.jl")
include("macros.jl")

# include("utils.jl")

end  # module Abstraction
