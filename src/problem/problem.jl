module Problem

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import MathematicalSystems as MS
import LazySets

using ..Utils
const UT = Utils

using ..System
const ST = System

include("horizon.jl")
include("problem_interface.jl")
include("control_problems.jl")
include("abstraction_problems.jl")
include("specifications.jl")
include("recipes.jl")

export OptimalControlProblem
export SafetyProblem
export ReachAndStayProblem
export CoSafeLTLProblem
export Infinity
export AbstractSpecification, StateSpec, TimedSpec, HybridSpec, hybrid_reach_spec

end
