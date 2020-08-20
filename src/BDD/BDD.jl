module BDD

using CUDD

# Helper functions

# Inline and remove assert later
function cube(manager::Ptr{CUDD.DdManager}, vars::Vector{Ptr{CUDD.DdNode}}, values::Vector{Cint})
    @assert length(vars) == length(values)
    return CUDD.Cudd_bddComputeCube(manager, vars, values, length(vars))
end

@inline _bit(x::T) where T<:Integer = iszero(x & one(T)) ? zero(Cint) : one(Cint)

function phase_rem(x::Int, n::Int)
    phase = zeros(Cint, n)
    for i in 1:n
        if !iszero(x & 1)
            phase[i] = one(Cint)
        end
        x >>= 1
    end
    return phase, x
end

include("BitSet.jl")
include("inttupleset.jl")

end # module
