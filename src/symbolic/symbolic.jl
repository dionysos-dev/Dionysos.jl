module Symbolic

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import ProgressMeter
import LazySets

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
include("grid_based_symbolic_model/local_symbolic_model.jl")
include("grid_based_symbolic_model/symbolic_model_list.jl")

# Abstraction computation: backend-agnostic orchestration + approximation kernels.
include("grid_based_symbolic_model/abstraction.jl")
include("grid_based_symbolic_model/transition_kernels.jl")

# Execution backends (how the transition relation is computed).
include("grid_based_symbolic_model/backends/backend.jl")
include("grid_based_symbolic_model/backends/threaded.jl")
include("grid_based_symbolic_model/backends/julia_distributed.jl")
include("grid_based_symbolic_model/backends/slurm_array.jl")

include("symbolic_hybrid_model/symbolic_hybrid_model.jl")

end  # module Symbolic
