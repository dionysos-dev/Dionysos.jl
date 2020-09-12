"""
    mutable struct IntTupleSet <: AbstractSet{NTuple{N,<:Integer} where N}

Same as `Base.Set{NTuple{N,<:Integer} where N}` but with `CUDD`.
"""
mutable struct IntTupleSet{N,T<:Integer} <: AbstractSet{NTuple{N,T}}
    mng::Ptr{Manager}
    root::Ptr{Node}
    indices_::NTuple{N,Vector{Cint}}
    auxindices_::NTuple{N,Vector{Cint}}
    phases1_::NTuple{N,Vector{Cint}}
    phases2_::NTuple{N,Vector{Cint}}
    auxphases_::NTuple{N,Vector{Cint}}
    values::Vector{Cint}
end

function IntTupleSet{N,T}() where {N,T}
    mng = CUDD.Cudd_Init()
    ARGS = ntuple(k -> ntuple(i -> Cint[], N), 5)
    root = _Zero(mng); _Ref(root)
    set = IntTupleSet{N,T}(mng, root, ARGS..., Cint[])
    finalizer(set) do set
        CUDD.Cudd_Quit(set.mng)
    end
    return set
end
IntTupleSet{N}() where N = IntTupleSet{N,Int}()

Base.eltype(::Type{IntTupleSet{N,T}}) where {N,T} = NTuple{N,T}
Base.empty(::IntTupleSet{N}, ::Type{T}=Int) where {N,T} = IntTupleSet{N,T}()
Base.emptymutable(::IntTupleSet{N}, ::Type{T}=Int) where {N,T} = IntTupleSet{N,T}()

# `x` will be used to denote a tuple of integers, and `e` to denote its elements

function _phases1!(set::IntTupleSet, x)
    for i in eachindex(set.indices_)
        indices = set.indices_[i]
        auxindices = empty!(set.auxindices_[i])
        phases = set.phases1_[i]
        auxphases = empty!(set.auxphases_[i])
        _compute_phases!(set.mng, phases, indices, auxphases, auxindices, x[i])
    end
end

# Returns the first `i` for which there is not enough bits to represent x[i];
# Returns 0 if there is no such `i`.
function _phases1_trunc!(set::IntTupleSet, x)
    for i in eachindex(set.indices_)
        !_compute_phases_trunc!(set.phases1_[i], x[i]) && return i
    end
    return 0
end

# Can we improve this?
@inline _incr_(e::T, i, j) where T = i == j ? zero(T) : (i == j + 1 ? e + one(T) : e)
function _increment(x::NTuple{N}, j) where N
    return ntuple(i -> _incr_(x[i], i, j), Val(N))
end

Base.iterate(set::IntTupleSet{N,T}) where {N,T} = iterate(set, ntuple(i -> zero(T), Val(N)))
function Base.iterate(set::IntTupleSet{N}, state::NTuple{N}) where N
    I = _phases1_trunc!(set, state)
    I == N && return nothing
    I == 0 && _eval_phases(set) && return (state, _increment(state, 0))
    return iterate(set, _increment(state, I))
end
