# ----------------------------
# Lazy set combinators
# ----------------------------

# Deduplicating iterator wrapper: makes `Set(itr)` / `enum_states` on a union
# reliable without materializing the whole set up front.
struct UniqueStates{I}
    iter::I
end

Base.IteratorEltype(::Type{UniqueStates{I}}) where {I} = Base.HasEltype()
Base.eltype(::Type{UniqueStates{I}}) where {I} = eltype(I)

Base.IteratorSize(::Type{UniqueStates{I}}) where {I} = Base.SizeUnknown()

function Base.length(U::UniqueStates)
    # O(n) but makes Set(itr) work reliably
    seen = BitSet()
    n = 0
    for v in U
        if !(v in seen)
            push!(seen, v)
            n += 1
        end
    end
    return n
end

function Base.iterate(U::UniqueStates, st = (iterate(U.iter), BitSet()))
    (it, seen) = st
    it === nothing && return nothing
    (val, s2) = it
    while in(val, seen)
        it = iterate(U.iter, s2)
        it === nothing && return nothing
        (val, s2) = it
    end
    push!(seen, val)
    return (val, (iterate(U.iter, s2), seen))
end

unique_states(iter) = UniqueStates(iter)

# ----------------------------

"""
    UnionStateSet{N, S1, S2} <: AbstractStateSet{N}

Lazy union `A ∪ B` of two state sets.
"""
struct UnionStateSet{N, S1 <: AbstractStateSet{N}, S2 <: AbstractStateSet{N}} <:
       AbstractStateSet{N}
    A::S1
    B::S2
end

contains_state(S::UnionStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    contains_state(S.A, m, q) || contains_state(S.B, m, q)

enum_states(S::UnionStateSet{N}, m::AbstractMapping{N}) where {N} =
    unique_states(Iterators.flatten((enum_states(S.A, m), enum_states(S.B, m))))

# ----------------------------

"""
    SetMinusStateSet{N, S1, S2} <: AbstractStateSet{N}

Lazy difference `A ∖ B` of two state sets.
"""
struct SetMinusStateSet{N, S1 <: AbstractStateSet{N}, S2 <: AbstractStateSet{N}} <:
       AbstractStateSet{N}
    A::S1
    B::S2
end

contains_state(S::SetMinusStateSet{N}, m::AbstractMapping{N}, q::Int) where {N} =
    contains_state(S.A, m, q) && !contains_state(S.B, m, q)

enum_states(S::SetMinusStateSet{N}, m::AbstractMapping{N}) where {N} =
    (q for q in enum_states(S.A, m) if !contains_state(S.B, m, q))

# `enum_states` is a lazy filter (no `length`), so the generic
# `get_n_state = length(enum_states(...))` cannot apply here; count directly.
get_n_state(S::SetMinusStateSet{N}, m::AbstractMapping{N}) where {N} =
    count(q -> !contains_state(S.B, m, q), enum_states(S.A, m))
