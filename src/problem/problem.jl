module Problem

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import MathematicalSystems as MS

using ..Utils
const UT = Utils

using ..System
const ST = System

include("problems.jl")

end
