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

include("symbolic_model.jl")

include("grid_based_symbolic_model/grid_based_symbolic_model.jl")
include("grid_based_symbolic_model/symbolic_model_list.jl")

include("symbolic_hybrid_model/symbolic_hybrid_model.jl")

end  # module Symbolic
