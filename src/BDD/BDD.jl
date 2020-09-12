module BDD

using CUDD

# Helper functions

# _One(mng::Ptr{Manager}) = CUDD.Cudd_ReadOne(mng)
_Deref(mng::Ptr{Manager}, node::Ptr{Node}) = Cudd_RecursiveDeref(mng, node)
_Ref(node::Ptr{Node}) = CUDD.Cudd_Ref(node)
_Cube(mng::Ptr{Manager}, vars::Vector{Ptr{Node}}, phases::Vector{Cint}) =
    CUDD.Cudd_bddComputeCube(mng, vars, phases, length(vars))
_Cube(mng::Ptr{Manager}, phases::Vector{Cint}) =
    CUDD.Cudd_IndicesToCube(mng, phases, length(values))
_Eval(mng::Ptr{Manager}, f::Ptr{Node}, values::Vector{Cint}) =
    CUDD.Cudd_Eval(mng, f, values) === CUDD.Cudd_ReadOne(mng)

@inline _bit(x::T) where T<:Integer = iszero(x & one(T)) ? zero(Cint) : one(Cint)

include("intset.jl")
include("inttupleset.jl")

Base.IteratorSize(::S) where {S<:Union{IntSet,IntTupleSet}} = Base.SizeUnknown()
Base.isempty(set::S) where {S<:Union{IntSet,IntTupleSet}} = set.root === set._ZERO
Base.empty!(set::IntSet) = _empty!(set)
Base.empty!(set::IntTupleSet) = _empty!(set)
function _empty!(set)
    _Deref(set.mng, set.root)
    set.root = set._ZERO
    _Ref(set.root)
    return set
end

Base.push!(set::IntSet{T}, x::T) where T = _push!(set, x)
Base.push!(set::IntTupleSet{N,T}, x::NTuple{N,T}) where {N,T} = _push!(set, x)
function _push!(set, x)
    _phase!(set, x)
    c = _Cube(set.mng, set.vars_, set.z_)
    _Ref(c)
    tmp1 = CUDD.Cudd_bddAnd(set.mng, set.root, c)
    _Ref(tmp1)
    _Deref(set.mng, set.root)
    _Deref(set.mng, c)
    c = _Cube(set.mng, set.variables, set.phases_)
    _Ref(c)
    tmp2 = CUDD.Cudd_bddOr(set.mng, tmp1, c)
    _Ref(tmp2)
    _Deref(set.mng, tmp1)
    _Deref(set.mng, c)
    set.root = tmp2
    return set
end

Base.in(x::T, set::IntSet{T}) where {N,T<:Integer} = _in(x, set)
Base.in(x, ::IntSet) = false
Base.in(x::NTuple{N,T}, set::IntTupleSet{N,T}) where {N,T<:Integer} = _in(x, set)
Base.in(x, ::IntTupleSet) = false
_in(x, set) = _phase_truncated!(set, x) == 0 && _Eval(set.mng, set.root, set.phases_)

Base.delete!(set::IntSet, x) = _delete!(set, x)
Base.delete!(set::IntTupleSet, x) = _delete!(set, x)
function _delete!(set, x)
    # ∈ updates `set.phases_`
    x ∈ set || return set
    # Use Nand because `Cudd_Not()` seems not implemented in CUDD
    c = _Cube(set.mng, set.variables, set.phases_)
    _Ref(c)
    notc = CUDD.Cudd_bddXor(set.mng, c, set._ONE)
    _Ref(notc)
    _Deref(set.mng, c)
    tmp = CUDD.Cudd_bddAnd(set.mng, set.root, notc)
    _Ref(tmp)
    _Deref(set.mng, notc)
    _Deref(set.mng, set.root)
    set.root = tmp
    return set
end

end # module
