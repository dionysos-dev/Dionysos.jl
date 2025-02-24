module System

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import MathematicalSystems as MS

using Symbolics
using IntervalArithmetic
using IntervalLinearAlgebra
import Plots: plot!

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
