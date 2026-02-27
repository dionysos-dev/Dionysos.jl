module Mapping

using Base.Iterators
import StaticArrays: SVector, SMatrix
import RecipesBase: @recipe, @series
import Colors

using ..Utils
const UT = Utils


@enum INCL_MODE INNER OUTER CENTER
@inline function _invInclMode(mode::INCL_MODE)
    mode === INNER && return OUTER
    mode === OUTER && return INNER
    mode === CENTER && return OUTER
    return error("Invalid inclusion mode $mode")
end


# ----------------------------
# Mapping API
# ----------------------------

"""
AbstractMapping{N,T}:
- defines a universe of labels (1..get_n_state)
- provides coord <-> state conversions
"""
abstract type AbstractMapping{N,T} end
get_n_state(m::AbstractMapping) = error("implement get_n_state")
enum_states(m::AbstractMapping) = 1:get_n_state(m)
get_state_by_coord(m::AbstractMapping, x) = error("implement get_state_by_coord")
get_coord_by_state(m::AbstractMapping, q::Int) = error("implement get_coord_by_state")
get_elem_by_state(m::AbstractMapping, q::Int) = error("implement get_elem_by_state")
get_elem_by_coord(m::AbstractMapping, x) = error("implement get_elem_by_coord")
enum_elems(m::AbstractMapping) = (get_elem_by_state(m, q) for q in enum_states(m))

# Dimension of the concrete space associated with the mapping
get_dim(::AbstractMapping{N,T}) where {N,T} = N
is_valid_state(m::AbstractMapping, q::Int) = (1 <= q <= get_n_state(m))
enum_coords(m::AbstractMapping) = (get_coord_by_state(m, q) for q in enum_states(m))

Base.empty!(m::AbstractMapping) =
    error("empty!(::$(typeof(m))) not implemented (mapping is read-only or not finite)")

get_states_from_set(::AbstractMapping, set, incl_mode::INCL_MODE) = error("implement get_states_from_set")

crop_to_mapping(m::AbstractMapping, states) = intersect(states, enum_states(m))
convert_to_list_mapping(m::AbstractMapping) =
    error("convert_to_list_mapping(::$(typeof(m))) not implemented")


@recipe function f(m::AbstractMapping{N,T}) where {N,T}
    S = MappingSet{N}()
    return ((S, m),) # delegates to the tuple recipe
end

include("list_mapping.jl")

include("grid_mapping/grid.jl")
include("grid_mapping/grid_mapping.jl")
include("grid_mapping/periodic_mapping.jl")

include("abstract_state_set.jl")

end
