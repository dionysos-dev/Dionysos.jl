module Symbolic

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import ProgressMeter

using LinearAlgebra, Colors
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
using ProgressMeter, IntervalArithmetic, LazySets
using JuMP
include("ellipsoidal_transitions.jl")

# nouveau code disparaitra à terme
include("HybridTimedAutomata/VectorContinuousSystem.jl")

include("HybridTimedAutomata/TimeSymbolicModel.jl")

include("HybridTimedAutomata/TimedHybridAutomata.jl")

end  # module Symbolic
