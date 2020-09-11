"""
    mutable struct IntSet <: AbstractSet{<:Integer}

Same as `Base.Set{<:Integer}` but with `CUDD`.
"""
mutable struct IntSet{T<:Integer} <: AbstractSet{T}
    mng::Ptr{Manager}
    variables::Vector{Ptr{Node}}
    root::Ptr{Node}
    phase_::Vector{Cint}
    vars_::Vector{Ptr{Node}}
    z_::Vector{Cint}
end

function IntSet{T}() where T
    mng = CUDD.Cudd_Init()
    variables = Ptr{Node}[]
    root = CUDD.Cudd_ReadLogicZero(mng)
    _Ref(root)
    phase_ = Cint[]
    z_ = Cint[]
    vars_ = Ptr{Node}[]
    set = IntSet{T}(mng, variables, root, phase_, vars_, z_)
    # finalizer(set) do
    #     CUDD.quit_cudd(set.mng)
    #     println("set desctructed")
    # end
    return set
end
IntSet() = IntSet{Int}()

Base.show(io::IO, set::IntSet) = Base.summary(io, set)
function Base.summary(io::IO, set::IntSet)
    nbits = length(set.variables)
    Base.showarg(io, set, true)
    print(io, " with ", nbits, " bit")
    isone(nbits) || print(io, "s")
end

Base.eltype(::Type{IntSet{T}}) where T = T
Base.empty(::IntSet, ::Type{T}=Int) where T = IntSet{T}()
Base.emptymutable(::IntSet, ::Type{T}=Int) where T = IntSet{T}()
Base.isempty(set::IntSet) = set.root === CUDD.Cudd_ReadLogicZero(set.mng)
function Base.empty!(set::IntSet)
    _Deref(set.mng, set.root)
    set.root = CUDD.Cudd_ReadLogicZero(set.mng)
    _Ref(set.root)
    return set
end
Base.IteratorSize(::IntSet) = Base.SizeUnknown()

function _phase!(set::IntSet, x)
    empty!(set.vars_)
    empty!(set.z_)
    for idx in eachindex(set.variables)
        set.phase_[idx] = _bit(x)
        x >>>= 1
    end
    while x > 0
        push!(set.phase_, _bit(x))
        # As the `mng` is only used by this struct, `bddIthVar` should be
        # the same as `bddNewVar`.
        newvar = CUDD.Cudd_bddIthVar(set.mng, length(set.variables))
        push!(set.variables, newvar)
        push!(set.vars_, newvar)
        push!(set.z_, zero(Cint))
        x >>>= 1
    end
end

function Base.push!(set::IntSet{T}, x::T) where T
    _phase!(set, x)
    c = _Cube(set.mng, set.vars_, set.z_)
    _Ref(c)
    tmp1 = CUDD.Cudd_bddAnd(set.mng, set.root, c)
    _Ref(tmp1)
    _Deref(set.mng, set.root)
    _Deref(set.mng, c)
    c = _Cube(set.mng, set.variables, set.phase_)
    _Ref(c)
    tmp2 = CUDD.Cudd_bddOr(set.mng, tmp1, c)
    _Ref(tmp2)
    _Deref(set.mng, tmp1)
    _Deref(set.mng, c)
    set.root = tmp2
    return set
end

function Base.delete!(set::IntSet{T}, x::T) where {T<:Integer}
    # ∈ updates `set.phase_`
    x ∈ set || return set
    # Use Nand because `Cudd_Not()` seems not implemented in CUDD
    c = _Cube(set.mng, set.variables, set.phase_)
    _Ref(c)
    notc = CUDD.Cudd_bddXor(set.mng, CUDD.Cudd_ReadOne(set.mng), c)
    _Ref(notc)
    _Deref(set.mng, c)
    tmp = CUDD.Cudd_bddAnd(set.mng, set.root, notc)
    _Ref(tmp)
    _Deref(set.mng, notc)
    _Deref(set.mng, set.root)
    set.root = tmp
    return set
end
Base.delete!(set::IntSet, x) = set

# Return 0 if there are enough bits to represent x.
# Returning a bool would have been more straighforward but this is to be consistent
# with IntTupleSet.
function _phase_truncated!(set::IntSet, x)
    for idx in eachindex(set.variables)
        set.phase_[idx] = _bit(x)
        x >>>= 1
    end
    return iszero(x) ? 0 : 1
end

function Base.in(x::T, set::IntSet{T}) where {N,T<:Integer}
    return _phase_truncated!(set, x) == 0 && _in(set.mng, set.root, set.phase_)
end
Base.in(x, ::IntSet) = false

function Base.iterate(set::IntSet{T}, state::T=zero(T)) where T
    _phase_truncated!(set, state) > 0 && return nothing
    _in(set.mng, set.root, set.phase_) && return (state, state + 1)
    return iterate(set, state + 1)
end
