module Control

using ..Utils
UT = Utils

using ..Domain
const D = Domain

using ..System
SY = System

using ..Symbolic
SY = Symbolic

using StaticArrays, LinearAlgebra, Plots

include("optimal_control.jl")
include("controller.jl")
include("lazy_abstraction.jl")
include("trajectory.jl")
end  # module control
