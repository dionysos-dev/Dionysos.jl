module BDD

using CUDD

mutable struct TupleUIntSet{N} <: AbstractSet{NTuple{N,UInt}}
    manager::Ptr{CUDD.DdManager}
    variables::NTuple{N,Vector{Ptr{CUDD.DdNode}}}
    root::Ptr{CUDD.DdNode}
    bitrep_::Vector{Cint}
    values_::Vector{Cint}
    vars1_::Vector{Ptr{CUDD.DdNode}}
    vars2_::Vector{Ptr{CUDD.DdNode}}
    function TupleUIntSet(dim)
        manager = CUDD.initialize_cudd()
        variables = ntuple(i -> Ptr{CUDD.DdNode}[], dim)
        root = CUDD.Cudd_ReadLogicZero(manager)
        bitrep_ = Cint[]
        values_ = Cint[]
        vars1_ = Ptr{CUDD.DdNode}[]
        vars2_ = Ptr{CUDD.DdNode}[]
        return new{dim}(manager, variables, root, bitrep_, values_, vars1_, vars2_)
    end
end

@inline _cube(manager, vars, values) = CUDD.Cudd_bddComputeCube(manager, vars, values, length(values))

# Use "x" for tuple, and "e" for its elements
@inline _bit(e) = iszero(e & 1) ? zero(Cint) : one(Cint)

# Compute full (not truncated) bitrep of integer e, and add variables to ith
# components to match the bitsize of e
function _bitrep_elem!(set, e, i)
    for var in set.variables[i]
        push!(set.bitrep_, _bit(e))
        push!(set.vars2_, var)
        e >>= 1
    end
    while e > 0
        push!(set.bitrep_, _bit(e))
        var = CUDD.Cudd_bddNewVar(set.manager)
        push!(set.variables[i], var)
        push!(set.values_, zero(Cint))
        push!(set.vars1_, var)
        e >>= 1
    end
end

function _bitrep_full!(set, x)
    empty!(set.bitrep_)
    empty!(set.values_)
    empty!(set.vars1_)
    empty!(set.vars2_)
    for (i, e) in enumerate(x)
        _bitrep_elem!(set, e, i)
    end
end

function Base.push!(set::TupleUIntSet{N}, x::NTuple{N,UInt}) where N
    _bitrep_full!(set, x)
    set.root = CUDD.Cudd_bddAnd(set.manager, set.root, _cube(set.manager, set.vars1_, set.values_))
    set.root = CUDD.Cudd_bddOr(set.manager, set.root, _cube(set.manager, set.vars2_, set.bitrep_))
end

function _bitrep_elem_trunc!(set, e, i)
    for var in set.variables[i]
        push!(set.bitrep_, _bit(e))
        push!(set.vars2_, var)
        e >>= 1
    end
    return iszero(e)
end

function _bitrep_full_trunc!(set, x)
    empty!(set.bitrep_)
    empty!(set.vars2_)
    for (i, e) in enumerate(x)
        !_bitrep_elem_trunc!(set, e, i) && return i
    end
    return 0
end

function _in(bitrep, set)
    return CUDD.Cudd_Eval(set.manager, set.root, bitrep) === CUDD.Cudd_ReadOne(set.manager)
end

function Base.in(x::NTuple{N,UInt}, set::TupleUIntSet{N}) where N
    return iszero(_bitrep_full_trunc!(set, x)) && _in(set.bitrep_, set)
end

function Base.isempty(set::TupleUIntSet)
    return set.root === CUDD.Cudd_ReadLogicZero(set.manager)
end

function Base.empty!(set::TupleUIntSet)
    set.root = CUDD.Cudd_ReadLogicZero(set.manager)
end

@inline _incspecial(e, i, j) = i == j ? zero(UInt) : (i == j + 1 ? e + one(UInt) : e)

function _increment(x::NTuple{N,UInt}, j) where N
    return ntuple(i -> _incspecial(x[i], i, j), Val(N))
end

function Base.iterate(set::TupleUIntSet{N},
        state::NTuple{N,UInt}=ntuple(i -> zero(UInt), Val(N))) where N
    idx_trunc = _bitrep_full_trunc!(set, state)
    idx_trunc == N && return nothing
    iszero(idx_trunc) && _in(set.bitrep_, set) && return (state, _increment(state, 0))
    return iterate(set, _increment(state, idx_trunc))
end

end # module BDD
