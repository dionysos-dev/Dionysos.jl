"""
    mutable struct IntTupleSet <: AbstractSet{NTuple{N,<:Integer} where N}

Same as `Base.Set{NTuple{N,<:Integer} where N}` but with `CUDD`.
"""
mutable struct IntTupleSet{N,T<:Integer} <: AbstractSet{NTuple{N,T}}
    mng::Ptr{Manager}
    variables::Vector{Ptr{Node}}
    indexes::NTuple{N,Vector{UInt16}}
    root::Ptr{Node}
    phases_::Vector{Cint}
    vars_::Vector{Ptr{Node}}
    z_::Vector{Cint}
    _ONE::Ptr{Node}
    _ZERO::Ptr{Node}
end

function IntTupleSet{N,T}() where {N,T}
    mng = CUDD.Cudd_Init()
    variables = Ptr{Node}[]
    indexes = ntuple(i -> UInt16[], N)
    phases_ = Cint[]
    z_ = Cint[]
    vars_ = Ptr{Node}[]
    _ONE = CUDD.Cudd_ReadOne(mng)
    _ZERO = CUDD.Cudd_ReadLogicZero(mng)
    _Ref(_ZERO)
    root = _ZERO
    _Ref(root)
    set = IntTupleSet{N,T}(mng, variables, indexes, root, phases_, vars_, z_, _ONE, _ZERO)
    finalizer(set) do set
        CUDD.Cudd_Quit(set.mng)
    end
    return set
end
IntTupleSet{N}() where N = IntTupleSet{N,Int}()

Base.show(io::IO, set::IntTupleSet) = Base.summary(io, set)
function Base.summary(io::IO, set::IntTupleSet)
    nbits = length.(set.indexes)
    Base.showarg(io, set, true)
    print(io, " with ", Base.dims2string(nbits), " bit")
    isone(prod(nbits)) || print(io, "s")
end

Base.eltype(::Type{IntTupleSet{N,T}}) where {N,T} = NTuple{N,T}
Base.empty(::IntTupleSet{N}, ::Type{T}=Int) where {N,T} = IntTupleSet{N,T}()
Base.emptymutable(::IntTupleSet{N}, ::Type{T}=Int) where {N,T} = IntTupleSet{N,T}()

# `x` will be used to denote a tuple of integers, and `e` to denote its elements

function _phase_simple!(set, e, i)
    for idx in set.indexes[i]
        set.phases_[idx] = _bit(e)
        e >>>= 1
    end
    while e > 0
        push!(set.phases_, _bit(e))
        # As the `mng` is only used by this struct, `bddIthVar` should be
        # the same as `bddNewVar`.
        newidx = length(set.variables)
        newvar = CUDD.Cudd_bddIthVar(set.mng, Cint(newidx))
        push!(set.variables, newvar)
        push!(set.indexes[i], newidx + 1)
        push!(set.vars_, newvar)
        push!(set.z_, zero(Cint))
        e >>>= 1
    end
end

function _phase!(set::IntTupleSet, x)
    empty!(set.vars_)
    empty!(set.z_)
    for (i, e) in enumerate(x)
        _phase_simple!(set, e, i)
    end
end

function _phase_simple_truncated!(set, e, i)
    for idx in set.indexes[i]
        set.phases_[idx] = _bit(e)
        e >>>= 1
    end
    return iszero(e)
end

# Returns the first `i` for which there is not enough bits to represent x[i];
# Returns 0 if there is no such `i`.
function _phase_truncated!(set::IntTupleSet, x)
    for (i, e) in enumerate(x)
        _phase_simple_truncated!(set, e, i) || return i
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
    I = _phase_truncated!(set, state)
    I == N && return nothing
    I == 0 && _Eval(set.mng, set.root, set.phases_) && return (state, _increment(state, 0))
    return iterate(set, _increment(state, I))
end
