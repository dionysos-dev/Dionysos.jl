"""
    mutable struct IntSet <: AbstractSet{<:Integer}

Same as `Base.Set{<:Integer}` but with `CUDD`.
"""
mutable struct IntSet{T<:Integer} <: AbstractSet{T}
    manager::Ptr{CUDD.DdManager}
    variables::Vector{Ptr{CUDD.DdNode}}
    root::Ptr{CUDD.DdNode}
    phase_::Vector{Cint}
    vars_::Vector{Ptr{CUDD.DdNode}}
    z_::Vector{Cint}
end

function IntSet{T}() where T
    manager = CUDD.initialize_cudd()
    variables = Ptr{CUDD.DdNode}[]
    root = CUDD.Cudd_ReadLogicZero(manager)
    phase_ = Cint[]
    z_ = Cint[]
    vars_ = Ptr{CUDD.DdNode}[]
    set = IntSet{T}(manager, variables, root, phase_, vars_, z_)
    # finalizer(set) do
    #     CUDD.quit_cudd(set.manager)
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
Base.isempty(set::IntSet) = set.root === CUDD.Cudd_ReadLogicZero(set.manager)
function Base.empty!(set::IntSet)
    set.root = CUDD.Cudd_ReadLogicZero(set.manager)
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
        # As the `manager` is only used by this struct, `bddIthVar` should be
        # the same as `bddNewVar`.
        newvar = CUDD.Cudd_bddIthVar(set.manager, length(set.variables))
        push!(set.variables, newvar)
        push!(set.vars_, newvar)
        push!(set.z_, zero(Cint))
        x >>>= 1
    end
end

function Base.push!(set::IntSet{T}, x::T) where T
    _phase!(set, x)
    set.root = CUDD.Cudd_bddAnd(set.manager, set.root, cube(set.manager, set.vars_, set.z_))
    set.root = CUDD.Cudd_bddOr(set.manager, set.root, cube(set.manager, set.variables, set.phase_))
    return set
end

function Base.delete!(set::IntSet{T}, x::T) where {T<:Integer}
    # ∈ updates `set.phase_`
    x ∈ set || return set
    # Use Nand because `Cudd_Not()` seems not implemented in CUDD
    set.root = CUDD.Cudd_bddAnd(set.manager, set.root,
        CUDD.Cudd_bddNand(set.manager, CUDD.Cudd_ReadOne(set.manager),
            cube(set.manager, set.variables, set.phase_)))
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
    return _phase_truncated!(set, x) == 0 && _in(set.manager, set.root, set.phase_)
end
Base.in(x, ::IntSet) = false

function Base.iterate(set::IntSet{T}, state::T=zero(T)) where T
    _phase_truncated!(set, state) > 0 && return nothing
    _in(set.manager, set.root, set.phase_) && return (state, state + 1)
    return iterate(set, state + 1)
end
