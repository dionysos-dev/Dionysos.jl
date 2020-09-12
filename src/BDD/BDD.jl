module BDD

using CUDD

# Helper functions

_One(mng::Ptr{Manager}) = CUDD.Cudd_ReadOne(mng)
_Zero(mng::Ptr{Manager}) = CUDD.Cudd_ReadLogicZero(mng)
_Deref(mng::Ptr{Manager}, node::Ptr{Node}) = Cudd_RecursiveDeref(mng, node)
_Ref(node::Ptr{Node}) = CUDD.Cudd_Ref(node)
# Remember that indexing starts at zero in C
_IthVar(mng::Ptr{Manager}, idx::Cint) = CUDD.Cudd_bddIthVar(mng, idx)
_Size(mng::Ptr{Manager}) = Cint(CUDD.Cudd_ReadSize(mng))
_AND(mng::Ptr{Manager}, f::Ptr{Node}, g::Ptr{Node}) = CUDD.Cudd_bddAnd(mng, f, g)
_OR(mng::Ptr{Manager}, f::Ptr{Node}, g::Ptr{Node}) = CUDD.Cudd_bddOr(mng, f, g)
_XOR(mng::Ptr{Manager}, f::Ptr{Node}, g::Ptr{Node}) = CUDD.Cudd_bddXor(mng, f, g)
_NOT(mng::Ptr{Manager}, f::Ptr{Node}) = _XOR(mng, f, _One(mng))
_Eval(mng::Ptr{Manager}, f::Ptr{Node}, values::Vector{Cint}) =
    CUDD.Cudd_Eval(mng, f, values) === _One(mng)
# _Cube(mng::Ptr{Manager}, vars::Vector{Ptr{Node}}, phases::Vector{Cint}) =
#     (@assert length(vars) === length(phases);
#     CUDD.Cudd_bddComputeCube(mng, vars, phases, length(vars)))
# _Cube(mng::Ptr{Manager}, phases::Vector{Cint}) =
#     CUDD.Cudd_IndicesToCube(mng, phases, length(values))

function _Cube(mng::Ptr{Manager}, indices::Vector{Cint}, phases::Vector{Cint})
    @assert length(indices) === length(phases)
    n = length(indices)
    cube = _One(mng); _Ref(cube)
    for r in n:-1:1 # In CUDD it is implemented in backward order, so it mimiched...
        var = _IthVar(mng, indices[r])
        if phases[r] === Cint(1)
            tmp = _AND(mng, var, cube); _Ref(tmp)
        elseif phases[r] === Cint(0)
            not_var = _NOT(mng, var); _Ref(not_var)
            tmp = _AND(mng, not_var, cube); _Ref(tmp)
            _Deref(mng, not_var)
        else
            throw("Phases must be a Cint equal to 0 or 1")
        end
        _Deref(mng, cube)
        cube = tmp
    end
    _Deref(mng, cube)
    return cube
end

@inline _bit(x::T) where T<:Integer = iszero(x & T(1)) ? Cint(0) : Cint(1)

include("intset.jl")
include("inttupleset.jl")

Base.IteratorSize(::S) where {S<:Union{IntSet,IntTupleSet}} = Base.SizeUnknown()
Base.isempty(set::S) where {S<:Union{IntSet,IntTupleSet}} = set.root === _Zero(set.mng)
Base.empty!(set::IntSet) = _empty!(set)
Base.empty!(set::IntTupleSet) = _empty!(set)
function _empty!(set)
    _Deref(set.mng, set.root)
    set.root = _Zero(set.mng); _Ref(set.root)
    return set
end

_dims2string(nbits) = length(nbits) == 1 ? nbits[1] : Base.dims2string(nbits)
Base.show(io::IO, set::S) where {S<:Union{IntSet,IntTupleSet}} = Base.summary(io, set)
function Base.summary(io::IO, set::S) where {S<:Union{IntSet,IntTupleSet}}
    nbits = length.(set.indices_)
    Base.showarg(io, set, true)
    print(io, " with ", _dims2string(nbits), " bit")
    isone(prod(nbits)) || print(io, "s")
end

#---
function _compute_phases!(mng, phases, indices, auxphases, auxindices, e)
    for r in eachindex(indices)
        phases[r] = _bit(e)
        e >>>= 1
    end
    nvars = _Size(mng)
    while e > 0
        _IthVar(mng, nvars)
        push!(phases, _bit(e))
        push!(indices, nvars)
        push!(auxphases, Cint(0))
        push!(auxindices, nvars)
        nvars += Cint(1)
        e >>>= 1
    end
end

Base.push!(set::IntSet{T}, x::T) where T = _push!(set, x)
Base.push!(set::IntTupleSet{N,T}, x::NTuple{N,T}) where {N,T} = _push!(set, x)
function _push!(set, x)
    _phases1!(set, x)
    mng = set.mng
    indices_ = set.indices_
    phases_ = set.phases1_
    auxindices_ = set.auxindices_
    auxphases_ = set.auxphases_
    root = set.root # root is already referenced
    for i in eachindex(indices_)
        cube = _Cube(mng, auxindices_[i], auxphases_[i]); _Ref(cube)
        tmp = _AND(mng, root, cube); _Ref(tmp)
        _Deref(mng, cube); _Deref(mng, root)
        root = tmp
    end
    # At this point root is referenced
    elem = _One(mng); _Ref(elem)
    for i in eachindex(indices_)
        cube = _Cube(mng, indices_[i], phases_[i]); _Ref(cube)
        tmp = _AND(mng, elem, cube); _Ref(tmp)
        _Deref(mng, cube); _Deref(mng, elem)
        elem = tmp
    end
    # At this point elem is referenced
    tmp = _OR(mng, root, elem); _Ref(tmp)
    _Deref(mng, root); _Deref(mng, elem)
    set.root = tmp
    return set
end

#---
# Returns true if e can be represented with bits in phases
function _compute_phases_trunc!(phases, e)
    for r in eachindex(phases)
        phases[r] = _bit(e)
        e >>>= 1
    end
    return iszero(e)
end

Base.in(x::T, set::IntSet{T}) where {N,T<:Integer} = _in(x, set)
Base.in(x, ::IntSet) = false
Base.in(x::NTuple{N,T}, set::IntTupleSet{N,T}) where {N,T<:Integer} = _in(x, set)
Base.in(x, ::IntTupleSet) = false
function _fill_values(values, phases_, indices_)
    for i in eachindex(indices_)
        indices = indices_[i]
        phases = phases_[i]
        for r in eachindex(indices)
            values[indices[r]+1] = phases[r]
        end
    end
end
function _eval_phases(set)
    resize!(set.values, _Size(set.mng))
    _fill_values(set.values, set.phases1_, set.indices_)
    return _Eval(set.mng, set.root, set.values)
end
function _in(x, set)
    _phases1_trunc!(set, x) !== 0 && return false
    return _eval_phases(set)
end

Base.delete!(set::IntSet, x) = _delete!(set, x)
Base.delete!(set::IntTupleSet, x) = _delete!(set, x)
function _delete!(set, x)
    # ∈ updates `set.phases1`
    x ∈ set || return set
    mng = set.mng
    indices_ = set.indices_
    phases_ = set.phases1_
    elem = _One(mng); _Ref(elem)
    for i in eachindex(indices_)
        cube = _Cube(mng, indices_[i], phases_[i]); _Ref(cube)
        tmp = _AND(mng, elem, cube); _Ref(tmp)
        _Deref(mng, cube); _Deref(mng, elem)
        elem = tmp
    end
    # At this point elem is referenced
    not_elem = _NOT(mng, elem); _Ref(not_elem)
    _Deref(mng, elem)
    tmp = _AND(mng, set.root, not_elem); _Ref(tmp)
    _Deref(mng, set.root); _Deref(mng, not_elem)
    set.root = tmp
    return set
end

end # module
