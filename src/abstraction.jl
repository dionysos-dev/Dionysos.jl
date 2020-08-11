module Abstraction

using LinearAlgebra
using ProgressMeter
using StaticArrays

@enum INCL_MODE INNER OUTER

include("rectangle.jl")
include("polyhedron.jl")
include("grid.jl")
include("domain.jl")
include("controlsystem.jl")
include("automaton.jl")
include("symbolicmodel.jl")
include("controller.jl")

end  # module Abstraction
