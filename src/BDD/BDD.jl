module BDD

using CUDD

# Helper functions

CUDD.Cudd_Init() = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0)

# Inline and remove assert later
function cube(manager::Ptr{CUDD.DdManager}, vars::Vector{Ptr{CUDD.DdNode}}, values::Vector{Cint})
    @assert length(vars) == length(values)
    return CUDD.Cudd_bddComputeCube(manager, vars, values, length(vars))
end

function _in(manager, root, phase)
    return CUDD.Cudd_Eval(manager, root, phase) === CUDD.Cudd_ReadOne(manager)
end

@inline _bit(x::T) where T<:Integer = iszero(x & one(T)) ? zero(Cint) : one(Cint)

include("intset.jl")
include("inttupleset.jl")

end # module
