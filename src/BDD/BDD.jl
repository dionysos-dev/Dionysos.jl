module BDD

using CUDD

# Helper functions

function cube(manager::Ptr{CUDD.DdManager}, vars::Vector{Ptr{CUDD.DdNode}}, values::Vector{Cint})
    @assert length(vars) == length(values)
    return CUDD.Cudd_bddComputeCube(manager, vars, values, length(vars))
end

function phase_rem(variables::Vector{Ptr{CUDD.DdNode}}, x::Int)
    phase = zeros(Cint, length(variables))
    for i in eachindex(phase)
        if !iszero(x & 1)
            phase[i] = one(Cint)
        end
        x >>= 1
    end
    return phase, x
end

include("BitSet.jl")

end # module
