# "This type of mathematical system is necessary to perform the Cartesian product between the dynamical system and time. Eventually, this code should be integrated into MathematicalSystems.jl and should no longer be present in the Dionysos codebase, or at least not in the symbolic part."

"""
    VectorContinuousSystem(systems::Vector{<:AbstractContinuousSystem})

Container for a vector of mathematical systems, allowing joint manipulation as a single system.
"""
struct VectorContinuousSystem <: AbstractContinuousSystem
    systems::Vector{AbstractContinuousSystem}
end

# Not defined for ConstrainedBlackBoxControlContinuousSystem and not particularly useful at the moment
stateset(sys::VectorContinuousSystem) =
    tuple((MathematicalSystems.stateset(s) for s in sys.systems)...)
inputset(sys::VectorContinuousSystem) =
    tuple((MathematicalSystems.inputset(s) for s in sys.systems)...)
statedim(sys::VectorContinuousSystem) =
    sum(MathematicalSystems.statedim(s) for s in sys.systems)
inputdim(sys::VectorContinuousSystem) =
    sum(MathematicalSystems.inputdim(s) for s in sys.systems)
