module System

import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Colors
import MathematicalSystems as MS

using Symbolics
import IntervalArithmetic as IA
import IntervalLinearAlgebra
# Workaround for IntervalLinearAlgebra not having `Interval` bound
if !isdefined(IntervalLinearAlgebra, :Interval)
    # bind IntervalLinearAlgebra.Interval to IntervalArithmetic.Interval
    IntervalLinearAlgebra.Interval = IA.Interval
end
IL = IntervalLinearAlgebra

import LinearAlgebra as LA
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
