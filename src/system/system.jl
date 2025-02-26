module System

using StaticArrays
using Symbolics
using MathematicalSystems
using IntervalArithmetic
using IntervalLinearAlgebra
import RecipesBase

using ..Utils
UT = Utils

include("controlsystem.jl")
include("controller.jl")
include("trajectory.jl")
end
