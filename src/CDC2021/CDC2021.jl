module CDC2021

using LinearAlgebra
using ProgressMeter
using StaticArrays
using Base.Cartesian
using HybridSystems

using ..Abstraction
const AB = Abstraction

include("general_domain.jl")
include("utils.jl")
include("lazy_abstraction.jl")
include("alternating_simulation.jl")
end  # module CDC2021
