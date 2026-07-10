module Symbolic

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import ProgressMeter

using JuMP

import LinearAlgebra as LA
using HybridSystems
import MathematicalSystems as MS

using ..Utils
const UT = Utils

using ..System
const ST = System

using ..Problem
const PR = Problem

using ..Mapping
const MP = Mapping

# The finite automaton is the abstraction's transition graph.
include("automata/automaton.jl")
include("automata/sorted_automaton_list.jl")
include("automata/indexed_automaton_list.jl")
include("automata/fast_indexed_automaton_list.jl")
include("automata/closed_loop.jl")

include("metadata.jl")
include("symbolic_model.jl")

include("grid_based_symbolic_model/grid_based_symbolic_model.jl")

include("grid_based_symbolic_model/execution_backend.jl")
include("grid_based_symbolic_model/sequential_threaded_backend.jl")
include("grid_based_symbolic_model/julia_distributed_backend.jl")
include("grid_based_symbolic_model/slurm_array_backend.jl")

include("grid_based_symbolic_model/symbolic_model_list.jl")

include("symbolic_hybrid_model/symbolic_hybrid_model.jl")

end  # module Symbolic
