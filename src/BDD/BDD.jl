module BDD

using CUDD

# Helper functions

_Deref(mng::Ptr{Manager}, node::Ptr{Node}) = Cudd_RecursiveDeref(mng, node)
_Ref(node::Ptr{Node}) = CUDD.Cudd_Ref(node)
_Cube(mng::Ptr{Manager}, vars::Vector{Ptr{Node}}, phases::Vector{Cint}) =
    CUDD.Cudd_bddComputeCube(mng, vars, phases, length(vars))
_Cube(mng::Ptr{Manager}, phases::Vector{Cint}) =
    CUDD.Cudd_IndicesToCube(mng, phases, length(values))
_in(mng::Ptr{Manager}, f::Ptr{Node}, values::Vector{Cint}) =
    CUDD.Cudd_Eval(mng, f, values) === CUDD.Cudd_ReadOne(mng)

@inline _bit(x::T) where T<:Integer = iszero(x & one(T)) ? zero(Cint) : one(Cint)

include("intset.jl")
include("inttupleset.jl")

end # module
