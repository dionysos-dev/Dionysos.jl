module Mapping

using Base.Iterators
import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series

using ..Domain
const DO = Domain

using ..Utils
const UT = Utils

include("mapping_continuous.jl")

end
