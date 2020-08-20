module Abstraction

using LinearAlgebra
using ProgressMeter
using StaticArrays
using Base.Cartesian

@enum INCL_MODE INNER OUTER

# See https://github.com/JuliaLang/julia/issues/37073
Base.hash(x::Tuple{}, h::UInt) = h + Base.tuplehash_seed
@generated function Base.hash(x::NTuple{N}, h::UInt) where N
    quote
        h += Base.tuplehash_seed
        @nexprs $N i -> h = hash(x[$N-i+1], h)
    end
end

include("rectangle.jl")
include("polyhedron.jl")
include("grid.jl")
include("domain.jl")
include("controlsystem.jl")
include("automaton.jl")
include("symbolicmodel.jl")
include("controller.jl")

end  # module Abstraction
