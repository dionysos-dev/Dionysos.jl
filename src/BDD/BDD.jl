module BDD

using CUDD

# Helper functions

function cube(manager::Ptr{CUDD.DdManager}, vars::Vector{Ptr{CUDD.DdNode}}, values::Vector{Cint})
    @assert length(vars) == length(values)
    return CUDD.Cudd_bddComputeCube(manager, vars, values, length(vars))
end

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

function phase_rem(x::UInt, n::UInt)
    phase = zeros(Cint, n)
    for i in 1:n
        if !iszero(x & 1)
            phase[i] = one(Cint)
        end
        x >>= 1
    end
    return phase, x
end

# include("BitSet.jl")
include("tupleset.jl")

end # module
