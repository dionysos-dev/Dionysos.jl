"""
    mode_matrices(system) -> Vector{Matrix{Float64}}

The linear map of every mode of a discrete switched system, in mode order.

`HybridSystems.discreteswitchedsystem` stores a mode either as a bare matrix or as a system
object carrying one in its `A` field; both are accepted so a switched system can be built
either way. Anything else is an error rather than a silent guess, since the matrices drive
the preimages a bisimulation is built from.

It lives here rather than in `System` because `PathCompleteFramework` needs it too, and `Utils`
is the only module both can see.
"""
function mode_matrices(system::HybridSystems.HybridSystem)
    maps = system.resetmaps
    A = Vector{Matrix{Float64}}(undef, length(maps))
    for (i, m) in enumerate(maps)
        if m isa AbstractMatrix
            A[i] = Array(m)
        elseif :A in fieldnames(typeof(m))
            A[i] = Array(getfield(m, :A))
        else
            error("Cannot extract a mode matrix from a resetmap of type $(typeof(m)).")
        end
    end
    return A
end
