module Relation

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

include("abstract_relation.jl")
# include("abstract_relation.jl")

end
