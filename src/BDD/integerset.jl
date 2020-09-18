# Manipulation of sets of integers

## IntSet
"""
    mutable struct IntSet <: AbstractSet{<:Integer}

Same as `Base.Set{<:Integer}` but with `CUDD`.
"""
mutable struct IntSet{T<:Integer} <: AbstractSet{T}
    mng::BDDManager
    vcs::Tuple{VariablesCluster}
    Rroot::Ref{Ptr{Node}}
end

function IntSet{T}(vc::VariablesCluster) where T
    mng = vc.mng
    root = _Zero(mng.dd); _Ref(root)
    Rroot = Ref(root)
    push!(vc.Rroots, Rroot)
    return IntSet{T}(mng, tuple(vc), Rroot)
end
IntSet(vc::VariablesCluster) = IntSet{Int}(vcs)

function IntSet{T}() where T
    mng = BDDManager()
    vc = VariablesCluster(mng)
    set = IntSet{T}(vc)
    finalizer(set -> finalize(set.mng), set)
    return set
end
IntSet() = IntSet{Int}()

Base.eltype(::Type{IntSet{T}}) where T = T
Base.empty(::IntSet, ::Type{T}=Int) where T = IntSet{T}()
Base.emptymutable(::IntSet, ::Type{T}=Int) where T = IntSet{T}()

_phases!(set::IntSet, x) = _phases_simple!(set.vcs[1], x)

# Returns 0 if there are enough bits to represent x.
# Returning a bool would have been more straighforward but this is to be consistent
# with IntTupleSet.
function _phases_trunc!(set::IntSet, x)
    not_trunc = _phases_trunc_simple!(set.vcs[1], x)
    return not_trunc ? 0 : 1
end



function Base.iterate(set::IntSet{T}, state::T=zero(T)) where T
    _phases_trunc!(set, state) > 0 && return nothing
    _eval_values_set(set) && return (state, state + 1)
    return iterate(set, state + 1)
end

## IntTupleSet
"""
    mutable struct IntTupleSet <: AbstractSet{NTuple{N,<:Integer} where N}

Same as `Base.Set{NTuple{N,<:Integer} where N}` but with `CUDD`.
"""
mutable struct IntTupleSet{N,T<:Integer} <: AbstractSet{NTuple{N,T}}
    mng::BDDManager
    vcs::NTuple{N,VariablesCluster}
    Rroot::Ref{Ptr{Node}}
end

function IntTupleSet{N,T}(vcs::NTuple{N,VariablesCluster}) where {N,T}
    mng = vcs[1].mng
    @assert all(vc -> vc.mng === mng, vcs)
    root = _Zero(mng.dd); _Ref(root)
    Rroot = Ref(root)
    for vc in vcs
        push!(vc.Rroots, Rroot)
    end
    return IntTupleSet{N,T}(mng, vcs, Rroot)
end
IntTupleSet{N}(vcs::NTuple{N,VariablesCluster}) where N = IntTupleSet{N,Int}(vcs)

function IntTupleSet{N,T}() where {N,T}
    mng = BDDManager()
    vcs = ntuple(i -> VariablesCluster(mng), Val(N))
    set = IntTupleSet{N,T}(vcs)
    finalizer(set -> finalize(set.mng), set)
    return set
end
IntTupleSet{N}() where N = IntTupleSet{N,Int}()

Base.eltype(::Type{IntTupleSet{N,T}}) where {N,T} = NTuple{N,T}
Base.empty(::IntTupleSet{N}, ::Type{T}=Int) where {N,T} = IntTupleSet{N,T}()
Base.emptymutable(::IntTupleSet{N}, ::Type{T}=Int) where {N,T} = IntTupleSet{N,T}()

function _phases!(set::IntTupleSet, x)
    for (vc, e) in zip(set.vcs, x)
        _phases_simple!(vc, e)
    end
end

# Returns the first i for which there is not enough bits to represent x[i];
# Returns 0 if there is no such i.
function _phases_trunc!(set::IntTupleSet, x)
    i = 1
    for (vc, e) in zip(set.vcs, x)
        !_phases_trunc_simple!(vc, e) && return i
        i += 1
    end
    return 0
end

# TODO: Can we improve this?
@inline _incr_(e::T, i, j) where T = i == j ? zero(T) : (i == j + 1 ? e + one(T) : e)
function _increment(x::NTuple{N}, j) where N
    return ntuple(i -> _incr_(x[i], i, j), Val(N))
end

Base.iterate(set::IntTupleSet{N,T}) where {N,T} = iterate(set, ntuple(i -> zero(T), Val(N)))
function Base.iterate(set::IntTupleSet{N}, state::NTuple{N}) where N
    I = _phases_trunc!(set, state)
    I === N && return nothing
    I === 0 && _eval_values_set(set) && return (state, _increment(state, 0))
    return iterate(set, _increment(state, I))
end

## High-level functions
Base.IteratorSize(::S) where {S<:Union{IntSet,IntTupleSet}} = Base.SizeUnknown()
Base.isempty(set::S) where {S<:Union{IntSet,IntTupleSet}} = set.Rroot[] === _Zero(set.mng.dd)
Base.empty!(set::IntSet) = _empty!(set)
Base.empty!(set::IntTupleSet) = _empty!(set)
function _empty!(set)
    _Deref(set.mng.dd, set.Rroot[])
    set.Rroot[] = _Zero(set.mng.dd); _Ref(set.Rroot[])
    return set
end

_dims2string(nbits) = length(nbits) === 1 ? nbits[1] : Base.dims2string(nbits)
Base.show(io::IO, set::S) where {S<:Union{IntSet,IntTupleSet}} = Base.summary(io, set)
function Base.summary(io::IO, set::S) where {S<:Union{IntSet,IntTupleSet}}
    nbits = map(vc -> vc.nvars, set.vcs)
    Base.showarg(io, set, true)
    print(io, " with ", _dims2string(nbits), " bit")
    isone(prod(nbits)) || print(io, "s")
end

_bit(x::T) where T<:Integer = iszero(x & T(1)) ? Cint(0) : Cint(1)

# Adds new variables to vc and updates the phases of vc so that they contain the
# full bit representation of e
function _phases_simple!(vc, e)
    for r in 1:vc.nvars
        vc.phases[r] = _bit(e)
        e >>>= 1
    end
    while e > 0
        add_var!(vc)
        vc.phases[vc.nvars] = _bit(e)
        e >>>= 1
    end
end

Base.push!(set::IntSet{T}, x::T) where T = _push!(set, x)
Base.push!(set::IntTupleSet{N,T}, x::NTuple{N,T}) where {N,T} = _push!(set, x)
function _push!(set, x)
    _phases!(set, x)
    # Here the phases of vcs[i] contain the full bit rep of x[i]
    dd = set.mng.dd
    elem = _One(dd); _Ref(elem)
    for vc in set.vcs
        cube = _Cube(dd, vc.indices, vc.phases); _Ref(cube)
        tmp = _AND(dd, elem, cube); _Ref(tmp)
        _Deref(dd, cube); _Deref(dd, elem)
        elem = tmp
    end
    # At this point elem is referenced
    tmp = _OR(dd, set.Rroot[], elem); _Ref(tmp)
    _Deref(dd, set.Rroot[]); _Deref(dd, elem)
    set.Rroot[] = tmp
    return set
end

# Returns true iff e can be represented with bits in phases
function _phases_trunc_simple!(vc, e)
    for r in 1:vc.nvars
        vc.phases[r] = _bit(e)
        e >>>= 1
    end
    return iszero(e)
end

Base.in(x::T, set::IntSet{T}) where {N,T<:Integer} = _in(x, set)
Base.in(x, ::IntSet) = false
Base.in(x::NTuple{N,T}, set::IntTupleSet{N,T}) where {N,T<:Integer} = _in(x, set)
Base.in(x, ::IntTupleSet) = false
function _eval_values_set(set)
    for vc in set.vcs
        set_values!(vc)
    end
    return _Eval(set.mng.dd, set.Rroot[], set.mng.values)
end
function _in(x, set)
    _phases_trunc!(set, x) !== 0 && return false
    # _phases_trunc updates the phases of the vc is set.vcs
    return _eval_values_set(set)
end

Base.delete!(set::IntSet, x) = _delete!(set, x)
Base.delete!(set::IntTupleSet, x) = _delete!(set, x)
function _delete!(set, x)
    x ∈ set || return set
    # ∈ updates the phases of the vc is set.vcs
    dd = set.mng.dd
    elem = _One(dd); _Ref(elem)
    for vc in set.vcs
        cube = _Cube(dd, vc.indices, vc.phases); _Ref(cube)
        tmp = _AND(dd, elem, cube); _Ref(tmp)
        _Deref(dd, cube); _Deref(dd, elem)
        elem = tmp
    end
    # At this point elem is referenced
    not_elem = _NOT(dd, elem); _Ref(not_elem)
    _Deref(dd, elem)
    tmp = _AND(dd, set.Rroot[], not_elem); _Ref(tmp)
    _Deref(dd, set.Rroot[]); _Deref(dd, not_elem)
    set.Rroot[] = tmp
    return set
end
