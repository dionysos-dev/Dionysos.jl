module System

using StaticArrays
using Symbolics
using MathematicalSystems
using IntervalArithmetic
using IntervalLinearAlgebra
using Plots

using ..Utils
UT = Utils

include("controlsystem.jl")
include("controller.jl")
include("trajectory.jl")

include("approximation/approximation.jl")
end
