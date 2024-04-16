module Symbolic

using HybridSystems
using StaticArrays
using LinearAlgebra

using ..Utils
const UT = Utils

using ..Domain
const DO = Domain

using ..System
const ST = System

include("automaton.jl")
include("symbolicmodel.jl")
include("hierarchical_symbolic.jl")
include("ellipsoidal_transitions.jl")
include("lazy_symbolic.jl")
include("alternating_simulation.jl")
include("proba_automaton.jl")

end  # module Symbolic
