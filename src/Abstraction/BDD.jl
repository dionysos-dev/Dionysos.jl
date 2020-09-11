_One(mng::Ptr{Manager}) = CUDD.Cudd_ReadOne(mng)
_Deref(mng::Ptr{Manager}, node::Ptr{Node}) = Cudd_RecursiveDeref(mng, node)
_Ref(node::Ptr{Node}) = CUDD.Cudd_Ref(node)
_Cube(mng::Ptr{Manager}, vars::Vector{Ptr{Node}}, phases::Vector{Cint}) =
    CUDD.Cudd_bddComputeCube(mng, vars, phases, length(vars))
_Cube(mng::Ptr{Manager}, phases::Vector{Cint}) =
    CUDD.Cudd_IndicesToCube(mng, phases, length(values))
_Eval(mng::Ptr{Manager}, f::Ptr{Node}, values::Vector{Cint}) =
    CUDD.Cudd_Eval(mng, f, values) === _One(mng)

@inline _bit(x::T) where T<:Integer = iszero(x & one(T)) ? zero(Cint) : one(Cint)

struct TheIndexes{N,M}
    sources::NTuple{N,Vector{UInt16}}
    labels::NTuple{M,Vector{UInt16}}
    targets::NTuple{N,Vector{UInt16}}
end

function TheIndexes{N,M}() where {N,M}
    sources = ntuple(i -> UInt16[], N)
    labels = ntuple(i -> UInt16[], M)
    targets = ntuple(i -> UInt16[], N)
    return TheIndexes(source, labels, targets)
end

mutable struct BDDManager{N,M}
    mng::Ptr{Manager}
    variables::Vector{Ptr{Node}}
    indexes::TheIndexes{N,M}
    phase_::Vector{Cint}
    vars_::Vector{Ptr{Node}}
    z_::Vector{Cint}
    Rroots::Vector{Ref{Ptr{Node}}}
end

function BDDManager{N,M}() where {N,M}
    mng = CUDD.Cudd_Init()
    variables = Ptr{Node}[]
    indexes = TheIndexes{N,M}()
    phase_ = Cint[]
    z_ = Cint[]
    vars_ = Ptr{Node}[]
    Rroots = Ref{Ptr{Node}}[]
    BDD = BDDManager{N,M}(mng, variables, indexes, phase_, vars_, z_, Rroots)
    finalizer(BDD) do BDD
        CUDD.Cudd_Quit(BDD.mng)
    end
    return BDD
end

function _add_new_variables(mng, Rroots, vars, phases)
    c = _Cube(mng, vars, phases)
    _Ref(c)
    for Rroot in Rroots
        tmp = CUDD.Cudd_bddAnd(mng, Rroot[], c)
        _Ref(tmp)
        _Deref(mng, Rroot[])
        Rroot[] = tmp
    end
    _Deref(mng, c)
end
