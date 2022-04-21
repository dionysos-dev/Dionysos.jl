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

include("automaton.jl")
include("symbolicmodel.jl")
include("ellipsoidal_transitions.jl")

include("lazy_symbolic.jl")
include("alternating_simulation.jl")
end  # module Symbolic
