"""
    mutable struct IntSet <: AbstractSet{<:Integer}

Same as `Base.Set{<:Integer}` but with `CUDD`.
"""
mutable struct IntSet{T<:Integer} <: AbstractSet{T}
    mng::Ptr{Manager}
    variables::Vector{Ptr{Node}}
    root::Ptr{Node}
    phases_::Vector{Cint}
    vars_::Vector{Ptr{Node}}
    z_::Vector{Cint}
end

function IntSet{T}() where T
    mng = CUDD.Cudd_Init()
    variables = Ptr{Node}[]
    root = CUDD.Cudd_ReadLogicZero(mng)
    _Ref(root)
    phases_ = Cint[]
    z_ = Cint[]
    vars_ = Ptr{Node}[]
    set = IntSet{T}(mng, variables, root, phases_, vars_, z_)
    finalizer(set) do set
        CUDD.Cudd_Quit(set.mng)
    end
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

function _phase!(set::IntSet, x)
    empty!(set.vars_)
    empty!(set.z_)
    for idx in eachindex(set.variables)
        set.phases_[idx] = _bit(x)
        x >>>= 1
    end
    while x > 0
        push!(set.phases_, _bit(x))
        # As the `mng` is only used by this struct, `bddIthVar` should be
        # the same as `bddNewVar`.
        newvar = CUDD.Cudd_bddIthVar(set.mng, length(set.variables))
        push!(set.variables, newvar)
        push!(set.vars_, newvar)
        push!(set.z_, zero(Cint))
        x >>>= 1
    end
end

# Return 0 if there are enough bits to represent x.
# Returning a bool would have been more straighforward but this is to be consistent
# with IntTupleSet.
function _phase_truncated!(set::IntSet, x)
    for idx in eachindex(set.variables)
        set.phases_[idx] = _bit(x)
        x >>>= 1
    end
    return iszero(x) ? 0 : 1
end

function Base.iterate(set::IntSet{T}, state::T=zero(T)) where T
    _phase_truncated!(set, state) > 0 && return nothing
    _Eval(set.mng, set.root, set.phases_) && return (state, state + 1)
    return iterate(set, state + 1)
end
