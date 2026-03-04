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

using ..System
const ST = System

using ..Problem
const PR = Problem

using ..Mapping
const MP = Mapping

include("automaton/automaton.jl")
include("automaton/sorted_automaton_list.jl")
include("automaton/indexed_automaton_list.jl")

include("symbolic_model.jl")

include("grid_based_symbolic_model/grid_based_symbolic_model.jl")
include("grid_based_symbolic_model/symbolic_model_list.jl")

include("symbolic_hybrid_model/symbolic_hybrid_model.jl")

end  # module Symbolic
