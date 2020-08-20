"""
    mutable struct IntTupleSet <: AbstractSet{NTuple{N,Int} where N}

Same as `Base.Set{NTuple{N,Int} where N}` but with `CUDD`.
"""
mutable struct IntTupleSet{N,T<:Integer} <: AbstractSet{NTuple{N,T}}
    manager::Ptr{CUDD.DdManager}
    variables::Vector{Ptr{CUDD.DdNode}}
    indexes::NTuple{N,Vector{UInt16}}
    root::Ptr{CUDD.DdNode}
    phase_::Vector{Cint}
    vars_::Vector{Ptr{CUDD.DdNode}}
    z_::Vector{Cint}
end

function IntTupleSet{N,T}() where {N,T}
    manager = CUDD.initialize_cudd()
    variables = Ptr{CUDD.DdNode}[]
    indexes = ntuple(i -> UInt16[], N)
    root = CUDD.Cudd_ReadLogicZero(manager)
    phase_ = Cint[]
    z_ = Cint[]
    vars_ = Ptr{CUDD.DdNode}[]
    return IntTupleSet{N,T}(manager, variables, indexes, root, phase_, vars_, z_)
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
Base.isempty(set::IntTupleSet) = set.root === CUDD.Cudd_ReadLogicZero(set.manager)
function Base.empty!(set::IntTupleSet)
    set.root = CUDD.Cudd_ReadLogicZero(set.manager)
    return set
end
Base.IteratorSize(::IntTupleSet) = Base.SizeUnknown()

# `x` will be used to denote a tuple of integers, and `e` to denote its elements

function _phase_simple!(set, e, i)
    for idx in set.indexes[i]
        set.phase_[idx] = _bit(e)
        e >>= 1
    end
    while e > 0
        push!(set.phase_, _bit(e))
        # As the `manager` is only used by this struct, `bddIthVar` should be
        # the same as `bddNewVar`.
        newidx = length(set.variables)
        newvar = CUDD.Cudd_bddIthVar(set.manager, Cint(newidx))
        push!(set.variables, newvar)
        push!(set.indexes[i], newidx + 1)
        push!(set.vars_, newvar)
        push!(set.z_, zero(Cint))
        e >>= 1
    end
end

function _phase!(set, x)
    empty!(set.vars_)
    empty!(set.z_)
    for (i, e) in enumerate(x)
        _phase_simple!(set, e, i)
    end
end

function Base.push!(set::IntTupleSet{N,T}, x::NTuple{N,T}) where {N,T}
    _phase!(set, x)
    set.root = CUDD.Cudd_bddAnd(set.manager, set.root, cube(set.manager, set.vars_, set.z_))
    set.root = CUDD.Cudd_bddOr(set.manager, set.root, cube(set.manager, set.variables, set.phase_))
    return set
end

function Base.delete!(set::IntTupleSet{N,T}, x::NTuple{N,T}) where {N,T<:Integer}
    # ∈ updates `set.phase_`
    x ∈ set || return set
    # Use Nand because `Cudd_Not()` seems not implemented in CUDD
    set.root = CUDD.Cudd_bddAnd(set.manager, set.root,
        CUDD.Cudd_bddNand(set.manager, CUDD.Cudd_ReadOne(set.manager),
            cube(set.manager, set.variables, set.phase_)))
    return set
end
Base.delete!(set::IntTupleSet, x) = set

function _phase_simple_truncated!(set, e, i)
    for idx in set.indexes[i]
        set.phase_[idx] = _bit(e)
        e >>= 1
    end
    return iszero(e)
end

function _phase_truncated!(set, x)
    for (i, e) in enumerate(x)
        _phase_simple_truncated!(set, e, i) || return i
    end
    return 0
end

function _in(phase, set)
    return CUDD.Cudd_Eval(set.manager, set.root, phase) === CUDD.Cudd_ReadOne(set.manager)
end

function Base.in(x::NTuple{N,T}, set::IntTupleSet{N,T}) where {N,T}
    return iszero(_phase_truncated!(set, x)) && _in(set.phase_, set)
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
    iszero(I) && _in(set.phase_, set) && return (state, _increment(state, 0))
    return iterate(set, _increment(state, I))
end
