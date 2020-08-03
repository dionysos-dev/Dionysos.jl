module Abstraction

using LinearAlgebra
using ProgressMeter

@enum INCL_MODE INNER OUTER

abstract type GridSpace{N} end

struct CellRef
    value::UInt64
end
# Called by `Base.hashindex`.
# We don't implement `hash(::CellRef, ::UInt64)` as it is not called by `Dict`
# and that would require hashing again.
Base.hash(ref::CellRef) = ref.value
Base.isless(ref1::CellRef, ref2::CellRef) = isless(ref1.value, ref2.value)

struct GridSpaceHash{N} <: GridSpace{N}
    orig::NTuple{N, Float64}
    h::NTuple{N, Float64}
    elems::Dict{CellRef, NTuple{N, Int}}
    overflow_ref::CellRef
    overflow_pos::NTuple{N, Int}
end

abstract type SubSet{N} end

struct SubSetHash{N} <: SubSet{N}
    grid_space::GridSpaceHash{N}
    elems::Set{CellRef}
end

struct ControlSystem{N, Fsys<:Function, Fbound<:Function}
    tstep::Float64
    sys_noise::NTuple{N, Float64}
    meas_noise::NTuple{N, Float64}
    sys_map::Fsys
    bound_map::Fbound
end

abstract type SymbolicModel end

mutable struct SymbolicModelHash{N1, N2, N3} <: SymbolicModel
    X_grid::GridSpaceHash{N1}
    U_grid::GridSpaceHash{N2}
    Y_grid::GridSpaceHash{N3}
    elems::Vector{Tuple{CellRef, CellRef, CellRef}}
    issorted::Bool
    isunique::Bool
end

include("rectangle.jl")
include("gridspace.jl")
include("subset.jl")
include("controlsystem.jl")
include("symbolicmodel.jl")
include("macros.jl")

# include("utils.jl")

end  # module Abstraction
