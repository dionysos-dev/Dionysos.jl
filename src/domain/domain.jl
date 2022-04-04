include("../utils/utils.jl")

module Domain

# using LinearAlgebra
# using ProgressMeter
using StaticArrays
# using Base.Cartesian
# using HybridSystems
using Main.Utils
UT = Main.Utils

@enum INCL_MODE INNER OUTER

include("grid.jl")
include("domain_list.jl")
end  # module Domain