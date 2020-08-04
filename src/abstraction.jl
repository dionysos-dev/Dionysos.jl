module Abstraction

using LinearAlgebra
using ProgressMeter

@enum INCL_MODE INNER OUTER

include("rectangle.jl")
include("gridspace.jl")
include("subset.jl")
include("controlsystem.jl")
include("automaton.jl")
include("symbolicmodel.jl")
include("controller.jl")

end  # module Abstraction
