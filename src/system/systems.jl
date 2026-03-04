
# ------------------------------
# MathematicalSystems Extension
# ------------------------------

# Used to perform the Cartesian product between dynamical systems and time dynamic
struct VectorContinuousSystem{T <: MS.AbstractContinuousSystem} <: MS.AbstractContinuousSystem
    systems::Vector{T}
end

MS.stateset(sys::VectorContinuousSystem) =
    Tuple(MS.stateset(s) for s in sys.systems)

MS.inputset(sys::VectorContinuousSystem) =
    Tuple(MS.inputset(s) for s in sys.systems)

MS.statedim(sys::VectorContinuousSystem)::Int =
    sum(MS.statedim(s) for s in sys.systems; init = 0)

MS.inputdim(sys::VectorContinuousSystem)::Int =
    sum(MS.inputdim(s) for s in sys.systems; init = 0)


# ------------------------------
# HybridSystems Extension
# ------------------------------
