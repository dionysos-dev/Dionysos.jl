module Domain

using StaticArrays, Plots
using ..Utils
const UT = Utils

"""
    abstract type DomainType{N, T} end

General abstract type for spatial domains in `N` dimensions with elements of type `T`.
"""
abstract type DomainType{N, T} end  # General abstract domain type

include("grid_domain/grid.jl")
include("grid_domain/grid_domain.jl")
include("grid_domain/domain_list.jl")
include("grid_domain/periodic_domain.jl")
include("grid_domain/general_domain.jl")
include("grid_domain/nested_domain.jl")

include("custom_domain.jl")
include("continuous_domain.jl")

end  # module Domain
