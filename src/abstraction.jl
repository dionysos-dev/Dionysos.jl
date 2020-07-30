module Abstraction

using LinearAlgebra
using StaticArrays

@enum INCL_MODE INNER OUTER

abstract type GridSpace end

struct GridSpaceHash{N} <: GridSpace
    dim::Int
    orig::SVector{N,Float64}
    h::SVector{N,Float64}
    elems::Dict{UInt64,SVector{N,Int}}
    overflow_ref::UInt64
    overflow_pos::SVector{N,Int}
end

abstract type SubSpace end

struct SubSpaceHash{N} <: SubSpace
    grid_space::GridSpaceHash{N}
    elems::Set{UInt64}
end

struct ControlSystem{N}
    tstep::Float64
    sys_noise::SVector{N,Float64}
    meas_noise::SVector{N,Float64}
    sys_map::Function
    bound_map::Function
end

abstract type TransitionMap end

struct TransitionMapHash{NX, NU, NY} <: TransitionMap
    X_grid::GridSpaceHash{NX}
    U_grid::GridSpaceHash{NU}
    Y_grid::GridSpaceHash{NY}
    elems::Vector{Tuple{UInt64,UInt64,UInt64}}
    infos::Dict{String,Bool}
end

# struct SymbolicModel
#     X_grid::GridSpace
#     U_grid::GridSpace
#     cont_sys::ControlSystem
#     trans_map::TransitionMap
#     infos::Dict{String,Bool}
# end

include("rectangle.jl")
include("gridspace.jl")
include("subspace.jl")
include("controlsystem.jl")
include("transitionmap.jl")
include("macros.jl")


# include("controlsystem.jl")
# include("plotting.jl")
# include("transitioncontroller.jl")



#=
mutable struct ControlSystem
    dim::Int
    tstep::Float64
    nsub::Int
    F_sys::Function
    sys_noise::Vector{Float64}
    meas_noise::Vector{Float64}
    F_bound::Function
    sys_map!::Function
    function ControlSystem(tstep, nsub, F_sys, sys_noise, meas_noise, F_bound)
        @assert length(sys_noise) == length(meas_noise)
        dim = length(sys_noise)
        cont_sys = new(dim, tstep, F_sys, sys_noise, meas_noise, F_bound, () -> ())
        _set_sys_map!(cont_sys)
        return cont_sys
   end
end

mutable struct SymbolicModel
    X_space::GridSpace
    U_space::GridSpace
    bound_map::Vector{Vector{Float64}}
    cont_sys::ControlSystem
    edge_list::Vector{Tuple{Int, Int, Int}}
    function SymbolicModel(X_space, U_space, cont_sys)
        @assert X_space.dim == cont_sys.dim
        return new(X_space, U_space, cont_sys, Vector{Tuple{Int}}[])
    end
end

mutable struct SubSpace
    grid_space::GridSpace
    idx_list::Vector{Int}
    function SubSpace(grid_space)
        return new(grid_space, Int[])
    end
end

mutable struct Controller
    sym_model::SymbolicModel
    enabled_inputs::Vector{Vector{Int}}
    function Controller(sym_model)
        return new(sym_model, Vector{Int}[])
    end
end

include("gridspace.jl")
include("controlsystem.jl")
include("symbolicmodel.jl")
include("subspace.jl")
include("controller.jl")
include("plotting.jl")
include("utils.jl")
=#

end  # module Abstraction
