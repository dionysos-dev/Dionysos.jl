"""
    ExplicitIdSet{N} <: AbstractStateSet{N}

A state set stored as an explicit `BitSet` of integer labels.
"""
mutable struct ExplicitIdSet{N} <: AbstractStateSet{N}
    bits::BitSet
end
ExplicitIdSet{N}() where {N} = ExplicitIdSet{N}(BitSet())

Base.copy(S::ExplicitIdSet{N}) where {N} = ExplicitIdSet{N}(copy(S.bits))

contains_state(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} = in(q, S.bits)
enum_states(S::ExplicitIdSet{N}, m::AbstractMapping{N}) where {N} = S.bits
add_state!(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} = push!(S.bits, q)
remove_state!(S::ExplicitIdSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    delete!(S.bits, q)
empty_states!(S::ExplicitIdSet{N}) where {N} = empty!(S.bits)
