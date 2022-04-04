module Control

using ..Utils
UT = Utils

using ..Symbolic
SY = Symbolic

# include("optimal_control.jl")
# include("q_learning.jl")
include("controller.jl")
end  # module control