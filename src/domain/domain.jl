module Domain

using StaticArrays, Plots
using ..Utils
UT = Utils

@enum INCL_MODE INNER OUTER CENTER

include("grid.jl")
include("domain_list.jl")
include("custom_domain.jl")
include("general_domain.jl")
include("nested_domain.jl")
end  # module Domain
