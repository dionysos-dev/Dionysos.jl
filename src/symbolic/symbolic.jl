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

include("grid_based_symbolic_model/grid_based_symbolic_model.jl")
include("grid_based_symbolic_model/symbolic_model_list.jl")
include("grid_based_symbolic_model/lazy_symbolic_model_list.jl")
include("grid_based_symbolic_model/hierarchical_symbolic.jl")

include("ellipsoidal_transitions.jl")
include("alternating_simulation.jl")
include("proba_automaton.jl")

end  # module Symbolic
