module Mapping

using Base.Iterators
import StaticArrays: SVector, SMatrix
import LinearAlgebra as LA
import RecipesBase: @recipe, @series
import LazySets

using ..Utils
const UT = Utils

# Inclusion modes live in Utils (shared with Problem); aliased here so the
# established `MP.INNER` / `MP.OUTER` / `MP.CENTER` spelling keeps working.
const INCL_MODE = UT.INCL_MODE
const INNER = UT.INNER
const OUTER = UT.OUTER
const CENTER = UT.CENTER
const invert_incl_mode = UT.invert_incl_mode

# ----------------------------
# Mapping API
# ----------------------------

"""
    AbstractMapping{N, T}

The concrete ↔ abstract discretization: a bijection between a universe of integer
state labels `1:get_n_state(m)` and points/cells of an `N`-dimensional concrete
space (element type `T`).

# Extending

Implement `get_n_state`, `get_state_by_coord`, `get_states_by_coord`,
`get_coord_by_state`, `get_elem_by_state`, `get_elem_by_coord`,
`has_overlapping_cells`. `enum_states`, `enum_coords`, `enum_elems`,
`is_valid_state`, `get_dim` have generic defaults.
"""
abstract type AbstractMapping{N, T} end
get_n_state(m::AbstractMapping) = error("implement `get_n_state` for $(typeof(m))")
enum_states(m::AbstractMapping) = 1:get_n_state(m)
get_state_by_coord(m::AbstractMapping, x) =
    error("implement `get_state_by_coord` for $(typeof(m))")
get_states_by_coord(m::AbstractMapping, x) =
    error("implement `get_states_by_coord` for $(typeof(m))")
has_overlapping_cells(m::AbstractMapping) =
    error("implement `has_overlapping_cells` for $(typeof(m))")

get_coord_by_state(m::AbstractMapping, q::Int) =
    error("implement `get_coord_by_state` for $(typeof(m))")
get_elem_by_state(m::AbstractMapping, q::Int) =
    error("implement `get_elem_by_state` for $(typeof(m))")
get_elem_by_coord(m::AbstractMapping, x) =
    error("implement `get_elem_by_coord` for $(typeof(m))")
enum_elems(m::AbstractMapping) = (get_elem_by_state(m, q) for q in enum_states(m))

# Dimension of the concrete space associated with the mapping
get_dim(::AbstractMapping{N, T}) where {N, T} = N
is_valid_state(m::AbstractMapping, q::Int) = (1 <= q <= get_n_state(m))
enum_coords(m::AbstractMapping) = (get_coord_by_state(m, q) for q in enum_states(m))

is_periodic(m::AbstractMapping) = false

Base.empty!(m::AbstractMapping) =
    error("empty!(::$(typeof(m))) not implemented (mapping is read-only or not finite)")

crop_to_mapping(m::AbstractMapping, states) = intersect(states, enum_states(m))
convert_to_list_mapping(m::AbstractMapping) =
    error("convert_to_list_mapping(::$(typeof(m))) not implemented")

@recipe function f(m::AbstractMapping{N, T}) where {N, T}
    S = FullStateSet{N}()
    return ((S, m),) # delegates to the tuple recipe
end

include("list_mapping.jl")

include("grid_mapping/grid.jl")
include("grid_mapping/grid_mapping.jl")
include("grid_mapping/explicit_grid_mapping.jl")
include("grid_mapping/implicit_grid_mapping.jl")
include("grid_mapping/periodic_grid_mapping.jl")
include("grid_mapping/multi_level_mapping.jl")

include("state_sets/state_set.jl")
include("state_sets/explicit_id_set.jl")
include("state_sets/full_state_set.jl")
include("state_sets/implicit_state_set.jl")
include("state_sets/combinators.jl")
include("state_sets/mapped_state_set.jl")
include("state_sets/recipes.jl")

end
