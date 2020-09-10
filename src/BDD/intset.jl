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
    return IntSet{T}(manager, variables, root, phase_, vars_, z_)
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
        x >>= 1
    end
    while x > 0
        push!(set.phase_, _bit(x))
        # As the `manager` is only used by this struct, `bddIthVar` should be
        # the same as `bddNewVar`.
        newvar = CUDD.Cudd_bddIthVar(set.manager, length(set.variables))
        push!(set.variables, newvar)
        push!(set.vars_, newvar)
        push!(set.z_, zero(Cint))
        x >>= 1
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
        x >>= 1
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

#=
phase_rem(x, set::IntSet) = phase_rem(x, length(set.variables))

function Base.iterate(set::IntSet, state::Int=0)
    phase, x = phase_rem(state, set)
    if iszero(x)
        if _in(phase, set)
            return (state, state + 1)
        else
            return iterate(set, state + 1)
        end
    else
        # state would require more variables than what is in `set.variables`
        # so it is necessarily not in `set`.
        return nothing
    end
end

function phase!(set::IntSet, x::Int)
    phase, x = phase_rem(x, set)
    while x > 0
        push!(phase, iszero(x & 1) ? zero(Cint) : one(Cint))
        # As the `manager` is only used by this struct, `bddIthVar` should be
        # the same as `bddAddVar`.
        var = CUDD.Cudd_bddIthVar(set.manager, Cint(length(set.variables)))
        push!(set.variables, var)
        set.root = CUDD.Cudd_bddAnd(set.manager, set.root, cube(set.manager, [var], [zero(Cint)]))
        x >>= 1
    end
    return phase
end
minterm!(set::IntSet, x::Int) = cube(set.manager, set.variables, phase!(set, x))
function Base.push!(set::IntSet, x::Int)
    # `minterm!` modifies `set.root` so we need to call it before the `bddOr` line where we access `set.root`.
    minterm = minterm!(set, x)
    set.root = CUDD.Cudd_bddOr(set.manager, set.root, minterm)
    return set
end
function _in(phase::Vector{Cint}, set::IntSet)
    return CUDD.Cudd_Eval(set.manager, set.root, phase) === CUDD.Cudd_ReadOne(set.manager)
end
function Base.in(x::Int, set::IntSet)
    phase, x = phase_rem(x, set)
    return iszero(x) && _in(phase, set)
end
=#
