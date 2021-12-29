module Abstraction

using LinearAlgebra
using ProgressMeter
using StaticArrays
using Base.Cartesian
using HybridSystems

@enum INCL_MODE INNER OUTER

include("sorted_vector_set.jl")
include("rectangle.jl")
include("polyhedron.jl")
include("grid.jl")
include("domain.jl")
include("controlsystem.jl")
include("automaton.jl")
include("symbolicmodel.jl")
include("controller.jl")

include("PWA/pwa.jl")
end  # module Abstraction
