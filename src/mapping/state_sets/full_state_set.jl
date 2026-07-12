"""
    FullStateSet{N} <: AbstractStateSet{N}

The read-only "all states" set: it contains every valid label of whatever
[`AbstractMapping`](@ref) it is queried against (`1:get_n_state(m)`). Used as the
default state set of a symbolic model ("consider every cell of the grid").

!!! note
    Do not confuse `FullStateSet` (a *state set* covering the whole universe) with
    [`MappedStateSet`](@ref) (a *bundle* of an arbitrary state set together with
    its mapping).
"""
struct FullStateSet{N} <: AbstractStateSet{N} end

Base.copy(::FullStateSet{N}) where {N} = FullStateSet{N}()

contains_state(::FullStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    is_valid_state(m, q)
enum_states(::FullStateSet{N}, m::AbstractMapping{N}) where {N} = 1:get_n_state(m)
get_n_state(::FullStateSet{N}, m::AbstractMapping{N}) where {N} = get_n_state(m)

add_state!(::FullStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("FullStateSet is read-only")
remove_state!(::FullStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    error("FullStateSet is read-only")
empty_states!(::FullStateSet{N}) where {N} = error("FullStateSet is read-only")
