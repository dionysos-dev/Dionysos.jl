module Symbolic

using HybridSystems
using StaticArrays
using LinearAlgebra

using ..Utils
UT = Utils

using ..Domain
DO = Domain

using ..System
ST = System

include("automaton.jl")
include("symbolicmodel.jl")
include("hierarchical_symbolic.jl")
include("ellipsoidal_transitions.jl")
include("lazy_symbolic.jl")
include("alternating_simulation.jl")
include("proba_automaton.jl")

end  # module Symbolic
