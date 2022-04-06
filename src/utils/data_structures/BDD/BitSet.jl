"""
    mutable struct BitSet <: AbstractSet{Int}

Same as `Base.BitSet` but with `CUDD`.
"""
mutable struct BitSet <: AbstractSet{Int}
    manager::Ptr{CUDD.DdManager}
    variables::Vector{Ptr{CUDD.DdNode}}
    root::Ptr{CUDD.DdNode}
    function BitSet()
        manager = CUDD.initialize_cudd()
        return new(manager, Ptr{CUDD.DdNode}[], CUDD.Cudd_ReadLogicZero(manager))
    end
end

Base.show(io::IO, set::BitSet) = Base.summary(io, set)
function Base.summary(io::IO, set::BitSet)
    n = length(set.variables)
    Base.showarg(io, set, true)
    print(io, " with ", n, " bit")
    if !isone(n)
        print(io, "s")
    end
end

Base.eltype(::Type{BitSet}) = Int
Base.empty(::BitSet, ::Type{Int}=Int) = BitSet()
Base.emptymutable(::BitSet, ::Type{Int}=Int) = BitSet()
Base.isempty(set::BitSet) = set.root === CUDD.Cudd_ReadLogicZero(set.manager)
function Base.empty!(set::BitSet)
    set.root = CUDD.Cudd_ReadLogicZero(set.manager)
    return set
end
Base.IteratorSize(::BitSet) = Base.SizeUnknown()

phase_rem(x, set::BitSet) = phase_rem(x, length(set.variables))

function Base.iterate(set::BitSet, state::Int=0)
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

function phase!(set::BitSet, x::Int)
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
minterm!(set::BitSet, x::Int) = cube(set.manager, set.variables, phase!(set, x))
function Base.push!(set::BitSet, x::Int)
    # `minterm!` modifies `set.root` so we need to call it before the `bddOr` line where we access `set.root`.
    minterm = minterm!(set, x)
    set.root = CUDD.Cudd_bddOr(set.manager, set.root, minterm)
    return set
end
function _in(phase::Vector{Cint}, set::BitSet)
    return CUDD.Cudd_Eval(set.manager, set.root, phase) === CUDD.Cudd_ReadOne(set.manager)
end
function Base.in(x::Int, set::BitSet)
    phase, x = phase_rem(x, set)
    return iszero(x) && _in(phase, set)
end
