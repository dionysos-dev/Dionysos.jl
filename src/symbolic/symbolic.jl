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
include("proba_automaton.jl")
include("markov_chain.jl")
include("nested_symbolic.jl")
end  # module Symbolic
