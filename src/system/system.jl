module System

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Colors
import MathematicalSystems as MS

using Symbolics
using IntervalArithmetic
using IntervalLinearAlgebra

import JuMP: MOI

using ..Utils
UT = Utils

include("controlsystem.jl")

using Colors
import HybridSystems
include("controller.jl")
include("trajectory.jl")

include("approximation/approximation.jl")
end
