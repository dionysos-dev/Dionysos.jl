module Symbolic

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import ProgressMeter

import LinearAlgebra as LA
using Colors
using HybridSystems, MathematicalSystems

using Graphs, SimpleWeightedGraphs

using ..Utils
const UT = Utils

using ..Domain
const DO = Domain

using ..System
const ST = System

include("automaton/automaton.jl")
include("automaton/sorted_automaton_list.jl")
include("automaton/indexed_automaton_list.jl")

include("grid_based_symbolic_model/grid_based_symbolic_model.jl")
include("grid_based_symbolic_model/symbolic_model_list.jl")
include("grid_based_symbolic_model/lazy_symbolic_model_list.jl")
include("grid_based_symbolic_model/hierarchical_symbolic.jl")

include("alternating_simulation.jl")
include("proba_automaton.jl")

using Polyhedra
import IntervalArithmetic as IA
using ProgressMeter, LazySets
using JuMP
include("ellipsoidal_transitions.jl")

include("timed_hybrid_automaton/vector_continuous_system.jl")

include("timed_hybrid_automaton/time_symbolic_model.jl")

include("timed_hybrid_automaton/symbolic_timed_hybrid_systems.jl")

end  # module Symbolic
