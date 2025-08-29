
# Container for a vector of mathematical systems, allowing joint manipulation as a single system.
# This type performs the Cartesian product between dynamical systems and time components.
struct VectorContinuousSystem{T <: AbstractContinuousSystem} <: AbstractContinuousSystem
    systems::Vector{T}
end

VectorContinuousSystem(systems::Vector{<:AbstractContinuousSystem}) =
    VectorContinuousSystem{AbstractContinuousSystem}(systems)

MathematicalSystems.stateset(sys::VectorContinuousSystem) =
    Tuple(MathematicalSystems.stateset(s) for s in sys.systems)

MathematicalSystems.inputset(sys::VectorContinuousSystem) =
    Tuple(MathematicalSystems.inputset(s) for s in sys.systems)

MathematicalSystems.statedim(sys::VectorContinuousSystem)::Int =
    sum(MathematicalSystems.statedim(s) for s in sys.systems; init = 0)

MathematicalSystems.inputdim(sys::VectorContinuousSystem)::Int =
    sum(MathematicalSystems.inputdim(s) for s in sys.systems; init = 0)
