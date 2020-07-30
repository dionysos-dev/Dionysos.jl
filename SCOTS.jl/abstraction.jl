module Abstraction

using LinearAlgebra
using StaticArrays
using PyPlot
using PyCall

@enum INCL_MODE INNER OUTER

abstract type GridSpace{N} end

struct GridSpaceHash{N} <: GridSpace{N}
    orig::NTuple{N, Float64}
    h::NTuple{N, Float64}
    elems::Dict{UInt64, NTuple{N, Int}}
    overflow_ref::UInt64
    overflow_pos::NTuple{N, Int}
end

abstract type SubSpace{N} end

mutable struct SubSpaceHash{N} <: SubSpace{N}
    grid_space::GridSpaceHash{N}
    elems::Vector{UInt64}
    issorted::Bool
    isunique::Bool
end

struct ControlSystem{N}
    tstep::Float64
    sys_noise::NTuple{N, Float64}
    meas_noise::NTuple{N, Float64}
    sys_map::Function
    bound_map::Function
end

abstract type SymbolicModel end

mutable struct SymbolicModelHash <: SymbolicModel
    X_grid::GridSpaceHash
    U_grid::GridSpaceHash
    Y_grid::GridSpaceHash
    elems::Vector{Tuple{UInt64, UInt64, UInt64}}
    issorted::Bool
    isunique::Bool
end

include("gridspace.jl")
include("subspace.jl")
include("controlsystem.jl")
include("symbolicmodel.jl")
include("macros.jl")

include("plotting.jl")
include("utils.jl")

end  # module Abstraction
