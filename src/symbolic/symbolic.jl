module Symbolic

using HybridSystems
using StaticArrays
using LinearAlgebra

using ..Utils
UT = Utils

using ..Domain
DO = Domain

using ..System
SY = System

# @enum INCL_MODE INNER OUTER

include("automaton.jl")
include("symbolicmodel.jl")
end  # module Symbolic